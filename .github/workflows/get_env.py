import platform
import sys
import os
import stat
from pathlib import Path


env_variable = "VIRTUAL_ENV"
env_variable_value = str(Path(sys.executable).absolute().parents[1])

if platform.system() == "Windows":
    ext = ".ps1"
    script = f'$Env:{env_variable}="{env_variable_value}"'
else:
    ext = ".sh"
    script = f'export {env_variable}="{env_variable_value}"'

script_file_path = Path(__file__).parent / f"set_virtual_env{ext}"

with open(script_file_path, "w") as f:
    f.write(script)


st = os.stat(script_file_path)
os.chmod(script_file_path, st.st_mode | stat.S_IEXEC)
