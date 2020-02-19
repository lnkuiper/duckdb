#include "duckdb/common/types/data_chunk.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/helper.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/common/serializer.hpp"
#include "duckdb/common/types/null_value.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"

using namespace duckdb;
using namespace std;

DataChunk::DataChunk() {
}

void DataChunk::InitializeEmpty(vector<TypeId> &types) {
	assert(types.size() > 0);
	for (index_t i = 0; i < types.size(); i++) {
		data.emplace_back(Vector(*this, types[i], nullptr));
	}
}

void DataChunk::Initialize(vector<TypeId> &types) {
	assert(types.size() > 0);
	InitializeEmpty(types);
	for (index_t i = 0; i < types.size(); i++) {
		data[i].Initialize();
	}
}

void DataChunk::Reset() {
	for (index_t i = 0; i < column_count(); i++) {
		data[i].Initialize();
	}
	SetCardinality(0);
}

void DataChunk::Destroy() {
	data.clear();
	SetCardinality(0);
}

Value DataChunk::GetValue(index_t col_idx, index_t index) const {
	assert(index < size());
	return data[col_idx].GetValue(sel_vector ? sel_vector[index] : index);
}

void DataChunk::SetValue(index_t col_idx, index_t index, Value val) {
	data[col_idx].SetValue(sel_vector ? sel_vector[index] : index, move(val));
}

void DataChunk::Reference(DataChunk &chunk) {
	assert(chunk.column_count() == column_count());
	SetCardinality(chunk);
	for (index_t i = 0; i < column_count(); i++) {
		data[i].Reference(chunk.data[i]);
	}
}

void DataChunk::Copy(DataChunk &other, index_t offset) {
	assert(column_count() == other.column_count());
	assert(other.size() == 0 && !other.sel_vector);

	for (index_t i = 0; i < column_count(); i++) {
		VectorOperations::Copy(data[i], other.data[i], offset);
	}
	other.SetCardinality(size() - offset);
}

void DataChunk::Append(DataChunk &other) {
	if (other.size() == 0) {
		return;
	}
	if (column_count() != other.column_count()) {
		throw OutOfRangeException("Column counts of appending chunk doesn't match!");
	}
	assert(!sel_vector);
	for (index_t i = 0; i < column_count(); i++) {
		VectorOperations::Append(other.data[i], data[i]);
	}
	SetCardinality(size() + other.size());
}

void DataChunk::ClearSelectionVector() {
	Normalify();
	if (!sel_vector) {
		return;
	}

	for (index_t i = 0; i < column_count(); i++) {
		data[i].ClearSelectionVector();
	}
	sel_vector = nullptr;
}

void DataChunk::Normalify() {
	for (index_t i = 0; i < column_count(); i++) {
		data[i].Normalify();
	}
}

vector<TypeId> DataChunk::GetTypes() {
	vector<TypeId> types;
	for (index_t i = 0; i < column_count(); i++) {
		types.push_back(data[i].type);
	}
	return types;
}

string DataChunk::ToString() const {
	string retval = "Chunk - [" + to_string(column_count()) + " Columns]\n";
	for (index_t i = 0; i < column_count(); i++) {
		retval += "- " + data[i].ToString() + "\n";
	}
	return retval;
}

void DataChunk::Serialize(Serializer &serializer) {
	// write the count
	serializer.Write<sel_t>(size());
	serializer.Write<index_t>(column_count());
	for (index_t i = 0; i < column_count(); i++) {
		// write the types
		serializer.Write<int>((int)data[i].type);
	}
	// write the data
	for (index_t i = 0; i < column_count(); i++) {
		auto type = data[i].type;
		if (TypeIsConstantSize(type)) {
			index_t write_size = GetTypeIdSize(type) * size();
			auto ptr = unique_ptr<data_t[]>(new data_t[write_size]);
			// constant size type: simple memcpy
			VectorOperations::CopyToStorage(data[i], ptr.get());
			serializer.WriteData(ptr.get(), write_size);
		} else {
			assert(type == TypeId::VARCHAR);
			// strings are inlined into the blob
			// we use null-padding to store them
			auto strings = (const char **)data[i].data;
			VectorOperations::Exec(sel_vector, size(), [&](index_t j, index_t k) {
				auto source = !data[i].nullmask[j] ? strings[j] : NullValue<const char *>();
				serializer.WriteString(source);
			});
		}
	}
}

void DataChunk::Deserialize(Deserializer &source) {
	auto rows = source.Read<sel_t>();
	index_t column_count = source.Read<index_t>();

	vector<TypeId> types;
	for (index_t i = 0; i < column_count; i++) {
		types.push_back((TypeId)source.Read<int>());
	}
	Initialize(types);
	// now load the column data
	SetCardinality(rows);
	for (index_t i = 0; i < column_count; i++) {
		auto type = data[i].type;
		if (TypeIsConstantSize(type)) {
			// constant size type: simple memcpy
			auto column_size = GetTypeIdSize(type) * rows;
			auto ptr = unique_ptr<data_t[]>(new data_t[column_size]);
			source.ReadData(ptr.get(), column_size);
			Vector v(*this, data[i].type, ptr.get());
			VectorOperations::ReadFromStorage(v, data[i]);
		} else {
			auto strings = (const char **)data[i].data;
			for (index_t j = 0; j < rows; j++) {
				// read the strings
				auto str = source.Read<string>();
				// now add the string to the StringHeap of the vector
				// and write the pointer into the vector
				if (IsNullValue<const char *>((const char *)str.c_str())) {
					strings[j] = nullptr;
					data[i].nullmask[j] = true;
				} else {
					strings[j] = data[i].AddString(str);
				}
			}
		}
	}
	Verify();
}

void DataChunk::MoveStringsToHeap(StringHeap &heap) {
	for (index_t c = 0; c < column_count(); c++) {
		if (data[c].type == TypeId::VARCHAR) {
			// move strings of this chunk to the specified heap
			auto source_strings = (const char **)data[c].GetData();
			auto old_buffer = move(data[c].buffer);
			if (data[c].vector_type == VectorType::CONSTANT_VECTOR) {
				data[c].buffer = VectorBuffer::CreateConstantVector(TypeId::VARCHAR);
				data[c].data = data[c].buffer->GetData();
				auto target_strings = (const char **)data[c].GetData();
				if (!data[c].nullmask[0]) {
					target_strings[0] = heap.AddString(source_strings[0]);
				}
			} else {
				data[c].buffer = VectorBuffer::CreateStandardVector(TypeId::VARCHAR);
				data[c].data = data[c].buffer->GetData();
				auto target_strings = (const char **)data[c].GetData();
				VectorOperations::ExecType<const char *>(data[c], [&](const char *str, index_t i, index_t k) {
					if (!data[c].nullmask[i]) {
						target_strings[i] = heap.AddString(source_strings[i]);
					}
				});
			}
		}
	}
}

void DataChunk::Hash(Vector &result) {
	assert(result.type == TypeId::HASH);
	VectorOperations::Hash(data[0], result);
	for (index_t i = 1; i < column_count(); i++) {
		VectorOperations::CombineHash(result, data[i]);
	}
}

void DataChunk::Verify() {
#ifdef DEBUG
	assert(size() <= STANDARD_VECTOR_SIZE);
	// verify that all vectors in this chunk have the chunk selection vector
	sel_t *v = sel_vector;
	for (index_t i = 0; i < column_count(); i++) {
		assert(data[i].sel_vector() == v);
		data[i].Verify();
	}
	// verify that all vectors in the chunk have the same count
	for (index_t i = 0; i < column_count(); i++) {
		assert(size() == data[i].size());
	}
#endif
}

void DataChunk::Print() {
	Printer::Print(ToString());
}
