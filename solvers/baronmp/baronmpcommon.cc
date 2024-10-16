#include "baronmpcommon.h"

#include <filesystem>
#include <cstdlib>
#include <vector>
#include <sstream>

#include <string>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include <cstring>
#include <csignal>


#ifdef _WIN32
#include <process.h>
#include <windows.h>
#else
#include <sys/wait.h>
#include <unistd.h>
#include <chrono>
#endif

namespace fs = std::filesystem;

namespace mp {

#ifdef WIN32
  constexpr const char* EXENAME = "baronin.exe"; 
#else
  constexpr const char* EXENAME = "baronin";
#endif
char SEP = fs::path::preferred_separator;

int BaronmpCommon::runBaron(const std::string& arg) {
  std::vector<std::string> args = { EXENAME, arg };
  return run(args);
}


std::string BaronmpCommon::make_cmdline(const std::vector<std::string>& args) {
  std::ostringstream cmdline;

  for (size_t i = 0; i < args.size(); ++i) {
    const std::string& arg = args[i];
    bool needs_quotes = arg.find_first_of(" \t") != std::string::npos;

    if (needs_quotes) {
      cmdline << '"';
    }

    for (char c : arg) {
      if (c == '"') {
        cmdline << "\\\"";  // Escape double quotes
      }
      else {
        cmdline << c;
      }
    }

    if (needs_quotes) {
      cmdline << '"';
    }

    if (i < args.size() - 1) {
      cmdline << ' ';
    }
  }

  return cmdline.str();
}


int BaronmpCommon::run(const std::vector<std::string>& args) {
  if (args.empty()) {
    fmt::print(stderr, "No command specified.\n");
    return -1;
  }

  std::string cmd_path = Utils::findFile(args[0]);
  if (cmd_path.empty()) {
    fmt::print(stderr, "Command not found: {}\n", args[0]);
    return -1;
  }

#ifdef _WIN32
  // Windows-specific implementation using CreateProcess
  std::string cmdline = make_cmdline(args);
  STARTUPINFO si;
  PROCESS_INFORMATION pi;
  ZeroMemory(&si, sizeof(si));
  si.cb = sizeof(si);
  ZeroMemory(&pi, sizeof(pi));

  // Create the process
  if (!CreateProcess(NULL, cmdline.data(), NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi)) {
    fmt::print(stderr, "CreateProcess failed ({})\n", GetLastError());
    return -1;
  }

  // Wait until the process exits
  WaitForSingleObject(pi.hProcess, INFINITE);

  // Get the exit code
  DWORD exit_code = 0;
  GetExitCodeProcess(pi.hProcess, &exit_code);

  // Close process and thread handles
  CloseHandle(pi.hProcess);
  CloseHandle(pi.hThread);

  return static_cast<int>(exit_code);

#else
  // Unix-like system implementation using fork and execve
  pid_t pid = fork();

  if (pid == -1) {
    fmt::print(stderr, "fork failed: {}\n", strerror(errno));
    return -1;
  }
  else if (pid == 0) {
    // Child process
    std::vector<char*> exec_args;
    for (const auto& arg : args) {
      exec_args.push_back(const_cast<char*>(arg.c_str()));
    }
    exec_args.push_back(nullptr);  // Null-terminate the argument list

    if (execve(cmd_path.c_str(), exec_args.data(), environ) == -1) {
      fmt::print(stderr, "execve failed: {}\n", strerror(errno));
      exit(EXIT_FAILURE);
      exit(EXIT_FAILURE);
    }
  }
  else {
    // Parent process: wait for the child to complete
    int status;
    if (waitpid(pid, &status, 0) == -1) {
      fmt::print(stderr, "waitpid failed: {}\n", strerror(errno));
      return -1;
    }

    if (WIFEXITED(status)) {
      return WEXITSTATUS(status);
    }
    else if (WIFSIGNALED(status)) {
      fmt::print(stderr, "Process terminated by signal {}\n", WTERMSIG(status));
      return -1;
    }
  }
#endif

  return 0;
}


/**
* Initialize baronDir and nlFilePath
*/
void BaronmpCommon::initDirectories(const std::string& stub, 
  const std::string& scratch, bool overwrite) {
  std::string td;


  try {
    initialDir = fs::current_path().string();
  }
  catch (const fs::filesystem_error& e) {
    fmt::print(stderr, "Error retrieving current directory: {}\n", e.what());
  }

  
  
  // Check if scratch directory is valid or create it if possible
  if (!scratch.empty()) {
    fs::path scratch_path(scratch);
    if (overwrite)
    {
      if (fs::exists(scratch_path))
        recrmdir(scratch_path.string());
    }
    if (!fs::exists(scratch_path)) {
      if (fs::create_directory(scratch_path)) {
        baronDir = scratch;
      }
    }
    else if (fs::is_directory(scratch_path)) {
      td = scratch;
    }
  }

  // Set up temporary directory
  if (baronDir.empty()) {
    if (td.empty())
    td = std::getenv("TEMP") ? std::getenv("TEMP") : "";
    if (td.empty() || !fs::exists(td) || !fs::is_directory(td)) {
#ifdef _WIN32
      fmt::print(stderr, "TEMP = \"{}\" is not a directory.\n", td);
      std::exit(1);
#else
      td = "/tmp";
#endif
    }
  // Build the temporary directory path
    fs::path tmp_dir = fs::path(td) / ("baron_tmp" + std::to_string(static_cast<unsigned long>(getpid())));
    if (fs::create_directory(tmp_dir)) {
      baronDir = tmp_dir.string();
    }
    else {
      fmt::print(stderr, "Could not create directory \"{}\".\n", tmp_dir.string());
      std::exit(1);
    }
  }

  // Construct nlfile path based on whether stub is absolute or relative
  fs::path nlfile_path;
  if (fs::path(stub).is_absolute()
#ifdef _WIN32
    || (stub.size() > 1 && stub[1] == ':')  // Also handle drive letter case on Windows
#endif
    ) {
    nlfile_path = fs::path(stub);
  }
  else {
    nlfile_path = fs::path(initialDir) / stub;
  }
  nlFilePath=nlfile_path.string();

  filePathBar = appendToDir(baronDir, FILENAME_BAR);
  filePathDic = appendToDir(baronDir, FILENAME_DIC);
  filePathAMPL= appendToDir(baronDir, FILENAME_AMPL);

  FILE_BAR = std::make_shared<fmt::BufferedFile>(filePathBar, "wb");
  std::filesystem::current_path(baronDir);

}

/**
* Recursively remove the directory and its contents
*/
int BaronmpCommon::recrmdir(const std::string& dname) {
  
  if (fs::is_directory(dname)) {
    
    std::error_code ec;
    fs::remove_all(dname, ec);
    if (ec) {
      fmt::print(stderr, "recrmdir: Failed to remove \"{}\": {}\n", dname, ec.message());
      return 1;
    }
    return 0;  // Success
  }
  else {
    fmt::print(stderr, "recrmdir: \"{}\" is not a directory.\n", dname);
    return 0;
  }
}

void BaronmpCommon::changeDirectory(const std::string& path) {
  try {
    fs::current_path(path);  // Set the current working directory
  }
  catch (const fs::filesystem_error& e) {
    fmt::print(stderr, "Error changing directory to '{}': {}\n", path, e.what());
  }
}


std::string BaronmpCommon::appendToDir(const std::string& dir, const std::string& file) {
  auto fullPath = fs::path(dir) / fs::path(file);
  return fullPath.string();
}

std::string ensureDoubleQuoted(const std::string& str) {
  if (str.size() >= 2 && str.front() == '"' && str.back() == '"') {
    return str;
  }
  std::string result = str;
  if (result.size() >= 2 && result.front() == '\'' && result.back() == '\'') {
    result = result.substr(1, result.size() - 2);
  }
  return fmt::format("\"{}\"", result);
}


void Options::initializeLPSolver() {
  const std::vector<std::string> avail = { "cbc", "cplex", "la04", "xpress" };
  if (std::find(avail.begin(), avail.end(), lpsolver) == avail.end())
    throw std::runtime_error(fmt::format("lpsolver {} not valid", lpsolver));

  lpsol = 8;
  if (lpsolver == "cplex")
  {
    lpsol = 3;
    if (cplexlibname.empty())
    {
      for (const auto& s : cplexlnames)
      {
        std::string res = Utils::findFile(s);
        if (!res.empty())
          cplexlibname = res;
        break;
      }
    }
  }
  else if (lpsolver == "xpress")
  {
    lpsol = 6;
    if (xpresslibname.empty())
    {
      for (const auto& s : xpresslnames)
      {
        std::string res = Utils::findFile(s);
        if (!res.empty())
          xpresslibname = res;
        break;
      }
    }
  }
}

void Options::appendInlineParams(fmt::MemoryWriter& m) const {
  if (inlineparams.size() > 0) {
    std::istringstream stream(inlineparams[0]);
    std::vector<std::string> tokens;
    std::string token;
    while (stream >> token) {
      tokens.push_back(token);
    }
    for (size_t i = 0; i < tokens.size(); i += 2)
      m << fmt::format("{}: {};\n", tokens[i], ensureDoubleQuoted(tokens[i + 1]));
  }
}

void BaronmpCommon::writeBaronOptions() {
  writeBaron("// AMPL/BARON default options\n");
  writeBaron("OPTIONS {\n");

  baronOptions().initializeLPSolver();
  writeBaron(baronOptions().serialize());

  writeBaron("}\n");
}
  // Write the globals data to a file
  bool GlobalsManager::writeGlobals(const std::string& filename, const  BaronGlobals& globals)   {
    try {
      std::ofstream output(filename);
      if (!output) {
        throw std::ios_base::failure("Failed to open file for writing.");
      }

      output << fmt::format("{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n",
        globals.baronDir, globals.initialDir,
        globals.lsolmsg, globals.lsolver, globals.lpsol_pvsub,  globals.nlfile,
        globals.nvars, globals.v_keepsol, globals.verbuf,
        globals.lpsol_dll);
      return true;
    }
    catch (const std::ios_base::failure& e) {
      fmt::print(stderr, "File write error: {}\n", e.what());
      return false;
    }
  }

  // Read the globals data from a file
  std::optional < GlobalsManager::BaronGlobals > GlobalsManager::readGlobals(const std::string& filename)  {
    BaronGlobals globals;
    try {
      std::ifstream input(filename);
      if (!input) {
        throw std::ios_base::failure("Failed to open file for reading.");
      }

      input >> globals.baronDir
        >> globals.initialDir
        >> globals.lsolmsg
        >> globals.lsolver
        >> globals.lpsol_pvsub
        >> globals.nlfile
        >> globals.nvars
        >> globals.v_keepsol
        >> globals.verbuf
        >> globals.lpsol_dll;

      if (input.fail()) {
        throw std::ios_base::failure("Error reading data from file.");
      }
      return globals;
    }
    catch (const std::ios_base::failure& e) {
      fmt::print(stderr, "File read error: {}\n", e.what());
      return std::nullopt;
    }
  }

  // Utility to display the contents of Globals
  void  GlobalsManager::displayGlobals(const  BaronGlobals& globals)  {
    fmt::print("Solver: {}\nLSolMsg: {}\nLPSubSolver: {}\nNumber of Vars: {}\nNLFile: {}\n"
      "Keep Solution: {}\nInitial Dir: {}\nBaron Dir: {}\nVerbBuf: {}\nLP DLL: {}\n",
      globals.lsolver, globals.lsolmsg, globals.lpsol_pvsub, globals.nvars,
      globals.nlfile, globals.v_keepsol, globals.initialDir, globals.baronDir,
      globals.verbuf, globals.lpsol_dll);
  }
  std::string Utils::getLibPath()
  {
    const std::string libName = "baron-lib";
    std::string libExtension;

#ifdef _WINDOWS
    libExtension = "dll";
#elif __APPLE__
    libExtension = "dylib";
#else
    libExtension = "so";
#endif

    return fmt::format("{}/{}.{}", getMyPath(true), libName, libExtension);
  }

std::string Utils::getExecutablePath() {
    #ifdef _WIN32
    // Windows
    char buffer[MAX_PATH];
    GetModuleFileName(nullptr, buffer, MAX_PATH);
    return std::string(buffer);

    #elif __APPLE__
    // macOS
    #include <mach-o/dyld.h>
    char buffer[PATH_MAX];
    uint32_t size = sizeof(buffer);
    if (_NSGetExecutablePath(buffer, &size) == 0) {
      return std::string(buffer);
    }
    return "";

    #elif __linux__
    // Linux
    #include <unistd.h>
    #include <limits.h>
    char buffer[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
    if (len != -1) {
      buffer[len] = '\0';
      return std::string(buffer);
    }
    return "";

    #else
    #error "Unsupported platform"
    #endif
}
std::string Utils::getMyPath(bool getDir, const std::string& appendPath) {
  // Get the full executable path
  fs::path executablePath = getExecutablePath();

  if (getDir) {
    auto pp = executablePath.parent_path();
    // Return the directory part of the path
    if (!appendPath.empty())
      pp /= appendPath;
    return pp.string();
  }
  else {
    // Return only the filename part of the path
    return executablePath.filename().string();
  }
}
void addFromString(const std::string& env_var,std::vector<std::string>& search_paths) {
  std::stringstream ss(env_var);
  std::string path;
  while (std::getline(ss, path, SEP)) {
    search_paths.emplace_back(path);
  }
}
std::string  Utils::findFile(const std::string& filename, bool isDLL) {
  
  // 0. If the filename is absolute or includes a directory separator, check it directly
  if (fs::path(filename).is_absolute() || filename.find_first_of("/\\") != std::string::npos) {
        return fs::exists(filename) ? filename : "";
  }

  // 1. Check current directory
  if (std::filesystem::exists(filename)) 
    return std::filesystem::absolute(filename).string();
  
  // 2. Check current exe directory
  auto current = getMyPath(true, filename);
  if (std::filesystem::exists(current))
    return current;

  std::vector<std::string> search_paths;
  if (const char* path_env = std::getenv("PATH"))
    addFromString(path_env, search_paths);
  if (const char* path_env = std::getenv("Path"))
    addFromString(path_env, search_paths);
  if (isDLL){
#ifdef __linux__
    if (const char* path_env = std::getenv("LD_LIBRARY_PATH"))
      addFromString(path_env, search_paths);

#endif
#ifdef __APPLE__
  if (const char* path_env = std::getenv("DYLD_LIBRARY_PATH"))
    addFromString(path_env, search_paths);
#endif
  }
  // 3. Search in specified directories
  for (const auto& path : search_paths) {
    std::filesystem::path lib_path = path;
    lib_path /= filename;
    if (std::filesystem::exists(lib_path)) {
      return std::filesystem::absolute(lib_path).string();
    }
  }
  // File not found
  return {};
}

TimFileData TimFileData::fromFile(const std::string& filePath) {
  TimFileData data;
  std::ifstream barTim(filePath, std::ios::in);
  std::string tmp;

  if (barTim.is_open()) {
    // Discard the first 5 values
    for (int i = 0; i < 5; ++i) {
      barTim >> tmp;
    }

    // Read printlb and printub values
    barTim >> tmp;
    data.printlb = Utils::readBaronValues(tmp).value();
    barTim >> tmp;
    data.printub = Utils::readBaronValues(tmp).value();

    // Read remaining status and timing values
    barTim >> data.barstatus >> data.modelstatus >> data.badmissingbounds;
    barTim >> data.itera >> data.nodeopt >> data.nodememmax;
    barTim >> data.totaltime;

    // Adjust itera if necessary
    if (data.itera < 0) {
      data.itera = 0;
    }
    barTim.close();
  }
  else {
    // Handle file error case
    data.Error_code = 3;  // Set error code if file can't be opened
  }

  return data;
}


/**
* Returns a tuple:
* Do read solution
* Baron terminated
* Start position
*/
std::tuple<bool, bool, int> findTermination(std::ifstream& barRes) {
  const std::string SOLSTRING = "***";
  const std::string BESTFOUND = "The best solution found";
  const std::string BETTERFOUND = "Better Solution Found";
  int pos;
  bool baron_solution_read=false, baron_terminated=false;
  int pos_solution=0;
  std::string line;
  //Find "***" or "best solution found" and get termination message
  while (!barRes.eof())
  {
    pos = barRes.tellg();
    getline(barRes, line);
    pos = line.find(SOLSTRING);
    if (pos != line.npos)
    {
      baron_solution_read = true;  //solution will be read
      pos_solution = pos;          // From here
      // TODO
      //Include termination message status to signature
      // interface_signature << line << endl;
      break;
    }

    //BARON run was terminated. Find for "best solution found"
    pos = line.find(BESTFOUND);
    if (pos != line.npos)
    {
      //Include termination message status to signature
      //interface_signature << line << endl;
      baron_solution_read = true; //solution will be read
      baron_terminated = true;
      //Finish current while loop
      break;
    }

    //BARON run was interrupted by user (barstatus = 6)
    //Find for "Better Solution Found", stores the last occurence
    pos = line.find(BETTERFOUND);
    if (pos != line.npos)
    {
      pos_solution = pos;
    }
  }
  return { baron_solution_read, baron_terminated, pos_solution };
}

// Static factory method to create an instance from file
ResFileData ResFileData::fromFile(const std::string& filePath, int nodeopt) {
  ResFileData data;
  std::ifstream barRes(filePath);
  std::string line;
  int barRes_position, pos, pos_betterfound;
  bool baron_solution_read;
  bool baron_terminated = false;
  // Only proceed if nodeopt is not -3
  if (nodeopt != -3 && barRes.is_open()) {
    // Read lines to parse necessary information
    while ((!data.numsol) && std::getline(barRes, line)) {
      // Parse for baron version, LP, NLP, and number of solutions
      data.parseBaronInfo(line);
    }
    std::tuple<bool, bool, int> res;

    // If numsol is not -1, skip until the termination message.
    // Otherwise, solutions will appear before the termination message.
    if (data.numsol != -1)
      res = findTermination(barRes);

     baron_terminated = std::get<1>(res);
    // Process all solutions if the file contains them
    
     
     if (baron_terminated)
       data.readSolution(barRes); // TODO
     else
     {
       bool more = false;
       do {
         more = data.readSolution(barRes);
       } while (more);
     }
  }
  else {
    // Handle case when file cannot be opened
    data.Error_code = 4;
  }
  barRes.close();
  return data;
}

std::optional<double> parseThirdItemAsDouble(const std::string& line) {
  std::istringstream iss(line);
  std::string item;
  if ((iss >> item) && (iss >> item) && (iss >> item)) {
    return Utils::readBaronValues(item);  // Successfully parsed, return the third item as double
  }
  return std::nullopt;
}

/**
* Scan the file for strings s1 (and s2 if not empty), and return 1 or 2 accordingly.
* Return 0 if not found.
* Also returns the line itself for further processing and the position of the seeked
* string.
* Note that the input stream position will be at the beginning of the following line,
* unless lookahead is specified 
*/
std::tuple<int, std::string, size_t> waitFor(std::ifstream& barRes, 
  const std::string& s1, const std::string& s2=std::string(), bool lookAhead=false) {
  size_t pos;
  std::string line;
  std::streampos spos;
  while (!barRes.eof())
  {
    if (lookAhead) spos = barRes.tellg();
    getline(barRes, line);
    
    pos = line.find(s1);
    if (pos != line.npos) {
      if (lookAhead) barRes.seekg(spos);
      return { 1, line, pos };
    }
    if (s2.empty()) continue;
    pos = line.find(s2);
    if (pos != line.npos)
    {
      if (lookAhead) barRes.seekg(spos);
      return { 2, line, pos };
    }
  }
  if (lookAhead) barRes.seekg(spos);
  return { 0, line, line.npos};

}

std::vector<double> readVector(std::ifstream& barRes, int size = 0) {
  std::string line = " ";
  std::vector<double> output;
  if (size) output.reserve(size);
  bool read = false;
  while (!barRes.eof() && !(line.length() == 0)) {
    getline(barRes, line);
    auto res = parseThirdItemAsDouble(line);
    if (res) {
      read = true;
      output.push_back(res.value());
    }
    else break;
  }
  if (!read)
    return std::vector<double>();
  return output;
}
bool ResFileData::readSolution(std::ifstream& barRes) {

  std::string line;
  std::string temp_string;
  bool found=false;
  
  auto res = waitFor(barRes, baron_objective);
  if (std::get<0>(res) == 0)
    return false;
  // line contains the current objective value
  line = std::get<1>(res);
  auto pos = std::get<2>(res);
  temp_string = line.substr(pos + strlen((char*)(baron_objective.c_str())));
  objective_values.push_back(Utils::readBaronValues(temp_string).value());
  // skip 2 lines
  getline(barRes, line); // Corresponding solution vector is:
  getline(barRes, line);  // Variable no.    Value
  // Read primal solution
  primal_solutions.push_back(readVector(barRes));

  res = waitFor(barRes, baron_objective, baron_dual, true);
  int resp = std::get<0>(res);
  switch (resp) {
  case 0: 
    return false;
  case 1:
    return true;
  case 2:
    getline(barRes, line); // Constraint no.  Price
    dual_solutions.push_back(readVector(barRes));
    // Or, if we read the duals, try to find the next header
    res = waitFor(barRes, baron_objective, "", true);
    return std::get<0>(res);
  }
  throw std::logic_error(fmt::format("Should not see this - waitFor = {}", resp));
}

// Helper function to parse baron-specific data in a line
void ResFileData::parseBaronInfo(const std::string& line) {
  size_t pos;
  if ((pos = line.find(baronVersion)) != std::string::npos) {
    baronVersion = line.substr(pos + baronVersion.size());
  }
  else if ((pos = line.find(baronLP)) != std::string::npos) {
    baronLP = line.substr(pos + baronLP.size());
  }
  else if ((pos = line.find(baronNLP)) != std::string::npos) {
    baronNLP = line.substr(pos + baronNLP.size());
  }
  else if ((pos = line.find(baronNumsol)) != std::string::npos) {
    numsol = std::stoi(line.substr(pos + baronNumsol.size()));
  }
}

} // namespace mp
