import json
import os
from pathlib import Path
import subprocess
import sys

TEST = 'A worker thread that cannot deschedule a blocked task tries again'
ROOT = Path(os.environ['GITHUB_WORKSPACE'])
OUT = ROOT / 'evidence'
OUT.mkdir(exist_ok=True)


def instrument(source):
    path = source / 'test/api/test_worker_exception.cpp'
    text = path.read_text()
    needle = '\tREQUIRE(blocked->shortest_gap >= RETRY_WAIT);'
    assert text.count(needle) == 1
    text = text.replace('#include <chrono>\n', '#include <chrono>\n#include <cstdio>\n')
    text = text.replace(needle, '\tstd::fprintf(stderr, "retry_wait_shortest_ns=%lld\\n",\n'
                        '\t             static_cast<long long>(std::chrono::duration_cast<std::chrono::nanoseconds>(blocked->shortest_gap).count()));\n' + needle)
    path.write_text(text)
    patch = subprocess.check_output(['git', 'diff', '--', str(path)], cwd=source, text=True)
    (OUT / (source.name + '-instrumentation.patch')).write_text(patch)


def run_test(source, suffix, args):
    binary = source / 'build/release/test/unittest.exe'
    result = subprocess.run([str(binary), *args], cwd=source, capture_output=True, text=True)
    output = result.stdout + result.stderr
    (OUT / (source.name + '-' + suffix + '.log')).write_text(output)
    gaps = [line for line in output.splitlines() if 'retry_wait_shortest_ns=' in line]
    print(source.name, suffix, 'returncode', result.returncode, 'gaps', gaps, flush=True)
    return {'case': suffix, 'returncode': result.returncode, 'gaps': gaps}


def build_and_test(name):
    source = ROOT / name
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=source, text=True).strip()
    print('SOURCE', name, revision, flush=True)
    instrument(source)
    env = os.environ.copy()
    env['EXTENSION_CONFIGS'] = (source / '.github/config/bundled_extensions.cmake').as_posix()
    subprocess.run(['make', 'release'], cwd=source, env=env, check=True)
    binary = source / 'build/release/test/unittest.exe'
    listed = subprocess.check_output([str(binary), '--list-tests'], cwd=source, text=True)
    (OUT / (name + '-all-tests.tsv')).write_text(listed)
    names = [line.split('\t')[0] for line in listed.splitlines() if line and line != 'name\tgroup']
    index = names.index(TEST)
    batch = names[(index // 10) * 10:(index // 10 + 1) * 10]
    batch_file = OUT / (name + '-original-batch.txt')
    batch_file.write_text('\n'.join('"' + test.replace('"', '\\"') + '"' for test in batch) + '\n')
    print('BATCH', name, index, batch, flush=True)
    results = []
    for i in range(10):
        results.append(run_test(source, 'isolated-' + str(i), [TEST]))
    for i in range(3):
        results.append(run_test(source, 'batch-' + str(i), ['-f', str(batch_file)]))
    summary = {'revision': revision, 'test_index': index, 'batch': batch, 'results': results}
    (OUT / (name + '-results.json')).write_text(json.dumps(summary, indent=2))


if __name__ == '__main__':
    build_and_test(sys.argv[1])
