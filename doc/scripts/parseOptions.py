import re
import os
import pandas as pd
import subprocess 
import sys
from collections import defaultdict

LOCAL_PATH =  os.path.dirname(os.path.abspath(__file__))
SOLVERS_PATH = os.path.join(LOCAL_PATH, "..", "..","..", "build", "vs64", "bin", "Debug")
solver_names = os.path.join(LOCAL_PATH,"solver_names.txt")

def is_executable(filepath):
    """Check if a file is an executable on both Windows and Unix-like systems."""
    if os.path.isfile(filepath):
        if sys.platform == "win32":
            return filepath.lower().endswith('.exe')
        else:
            return os.access(filepath, os.X_OK)
    return False

def find_executables(directory) -> dict:
    """Find all executables in the given directory."""
    executables = {}
    for root, _, files in os.walk(directory):
        for file in files:
            filepath = os.path.join(root, file)
            if is_executable(filepath):
                executables[os.path.splitext(os.path.basename(filepath))[0]]=filepath
    return executables

def run_executable(executable):
    """Run the executable with '-=' and capture the output."""
    try:
        result = subprocess.run([executable, "-="], capture_output=True, text=True)
        return result.stdout
    except subprocess.SubprocessError as e:
        print(f"Error running {executable}: {e}")
        return ""

 
 
def create_backend_mapping(solver_files: dict, backend_mapping_file: str):
    backend_names_mapping={}
    # Read backend mapping from the CSV file
    with open(backend_mapping_file, 'r') as file:
        for line in file:
            backend_name, url = line.strip().split(',')
            # Regular expression to match the word after 'solvers/' and before the next '/'
            match = re.search(r'solvers/(\w+)/', url)

            # Extract and print the matched word
            if match:
                friendly_name = match.group(1)
            backend_names_mapping[backend_name] = friendly_name

    # Filter solvers names based on the mapping file
    filtered_executables = { backend : solver_files[backend] for backend in backend_names_mapping.values() if backend in solver_files.keys()}
    return (backend_names_mapping, filtered_executables)


def parse_output(output):
    """Parse the output to extract options, aliases, and descriptions."""
    options = []
    option_pattern = re.compile(
        r'(\w+:\w+)\s+\(([^)]+)\)\s*\n\s*(.+?)(?=\n\s*\w+:\w+|\Z)',
        re.DOTALL
    )
    matches = option_pattern.finditer(output)
    for match in matches:
        main_name = match.group(1)
        aliases = match.group(2).split(", ")
        description = match.group(3).strip()

        if filter_option_name(main_name):
            options.append({
                'main_name': main_name,
                'aliases': aliases,
                'description': description
            })
    return options
    

def filter_option_name(name: str):
    """return true if the option has to be considered"""
    return not (name.startswith("acc:") or name.startswith("cvt:") or name.startswith("chk:") )
def group_common_options(executable_files):
    """Group options by main_name and list the executables that support each option."""
    common_options = defaultdict(list)
    solver_options = {}
    for name, exe in executable_files.items():
        output = run_executable(exe)
        options = parse_output(output)
        solver_options[name]=options
        for option in options:
            common_options[option['main_name']].append(name)
    
    return solver_options, common_options

def filter_common_options(options: dict, min_num: int):
    return {o: s for o,s in options.items() if len(s) == min_num}


executable_files = find_executables(SOLVERS_PATH)
(backends_map, executable_files) =create_backend_mapping(executable_files, solver_names)
all_solver_options, common_options = group_common_options(executable_files)


def find_not_supported(options: dict, solverlist: list) -> list:
    return {o: [s for s in solverlist if s not in supported] for (o, supported) in options.items() }


def merge_supported_not_supported(supported: dict, not_supported: dict) -> pd.DataFrame:
    df_supported = pd.DataFrame({
        'option': list(supported.keys()),
        'supported_by': list(supported.values())
    })
    df_not_supported = pd.DataFrame({
        'option': list(not_supported.keys()),
        'not_supported_by': list(not_supported.values())
    })
    df_merged = pd.merge(df_supported, df_not_supported, on='option', how='outer')
    df_merged['supported_by'] = df_merged['supported_by'].apply(lambda x: x if isinstance(x, list) else [])
    df_merged['not_supported_by'] = df_merged['not_supported_by'].apply(lambda x: x if isinstance(x, list) else [])
    return df_merged

print(all_solver_options)
for i in range(len(executable_files),1, -1):
    print(f"Supported by {i} solvers out of {len(executable_files)}")
    supported = filter_common_options(common_options, i)
    not_supported = find_not_supported(supported, executable_files.keys())
    df = merge_supported_not_supported(supported, not_supported)
    print(df)  
