import os
import sys
import subprocess

def local_branch_exists(branch):
    result = subprocess.run(
        ['git', 'branch', '--list', branch],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    return result.stdout.strip() != ""

if len(sys.argv) < 2:
    print("no new branch was given")
    raise Exception
    
remote_branch_name = sys.argv[1]
print("switch branch to",remote_branch_name)

if local_branch_exists(remote_branch_name):
    print("Branch exists locally: switching normally")
    subprocess.run(['git', 'checkout', remote_branch_name])
else:
    print("Branch does not exist locally: fetch from origin")
    subprocess.run(['git', 'remote', 'set-branches', 'origin', str(remote_branch_name)])
    subprocess.run(['git', 'fetch', '-v', '--depth=1'])
    subprocess.run(['git', 'checkout', '-b', str(remote_branch_name), 'origin/'+str(remote_branch_name)])

subprocess.run(['git', 'submodule', 'update', '--recursive', '--init'])
subprocess.run(['cmake', '.'])

if len(sys.argv) < 3:
    print("default to make install")
    doInstall = True
else:
    print("build CPlantBox locally")
    doInstall = False


if doInstall:
    subprocess.run(['make', 'install'])
else:
    subprocess.run(['make'])