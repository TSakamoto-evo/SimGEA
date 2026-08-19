#include "options.hpp"

std::string print_version(){
  std::string ver = "SimGEA neutral simulations v1.1.0";
  std::cout << "Version: " << ver << "\n";

  return(ver);
}

bool parse_int_c(const char* s, int& out){
  if(!s){
    return false;
  }

  errno = 0;
  char* end = nullptr;
  long v = std::strtol(s, &end, 10);

  if(errno == ERANGE || end == s || *end != '\0'){
    return false;
  }
  if(v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max()){
    return false;
  }

  out = static_cast<int>(v);
  return true;
}

bool parse_double_c(const char* s, double& out){
  if(!s){
    return false;
  }

  errno = 0;
  char* end = nullptr;
  double v = std::strtod(s, &end);

  if(errno == ERANGE || end == s || *end != '\0'){
    return false;
  }
  if(std::isnan(v) || !std::isfinite(v)){
    return false;
  }

  out = v;
  return true;
}

bool parse_uint_c(const char* s, unsigned int& out){
  if(!s){
    return false;
  }

  if(*s == '\0' || *s == '-'){
    return false;
  }

  errno = 0;
  char* end = nullptr;
  unsigned long v = std::strtoul(s, &end, 10);

  if(errno == ERANGE || end == s || *end != '\0'){
    return false;
  }
  if(v > std::numeric_limits<unsigned int>::max()){
    return false;
  }

  out = static_cast<unsigned int>(v);
  return true;
}


void print_help(const char* prog){
  print_version();
  std::cout << "\nUsage: " << prog << " [options] -m <migration_file> -d <baypass_file> -s <site_number> \n\n"

    << "Required:\n"
    << "  -m, --mig <STR>             Migration file (space-separated).\n"
    << "  -d, --data <STR>            Baypass style file (space-separated).\n"
    << "  -s, --site <INT>            Total number of variant sites per window\n\n"
    << "Optional (defaults in [brackets]):\n"
    << "  -o, --output <STR>          Output simulated file [simulated.txt]\n"
    << "      --mode <INT>            0: simulation mode, 1: simulation + GEA mode [0]\n"
    << "      --sim-mode <INT>        0: window-based, 1 bin-based (allowed in mode 0) [0]\n"
    << "  -b  --bin-file <STR>        File specifying bin breaks [Required for bin-based mode]\n"
    << "  -w, --win-width <DOUBLE>    Width of frequency window [0.025]\n"
    << "      --prerun <INT>          Number of trees in the prerun [10000]\n"
    << "      --mut-density <DOUBLE>  Expected number of mutations per window per tree [0.1]\n"
    << "      --seed <UNSIGNED INT>   Random seed [randomly determined]\n"
    << "  -l, --log-file <STR>        Log file [log.txt]\n"
    << "  -h, --help                  Show this help and exit\n\n"
    << "Options for +GEA mode (defaults in [brackets]):\n"
    << "  -e, --env <STR>             Baypass style environment file [Required]\n"
    << "      --pval <STR>            Output p-values file [p_vals.txt]\n"
    << "      --q-tail <DOUBLE>       The definition of the tail in the generalized Pareto fitting [0.995]\n"
    << std::endl;    
}

Options parse_args(int argc, char* argv[]){
  Options opt;

  if(argc <= 1){
    print_help(argv[0]);
    std::exit(1);
  }

  for(int i = 1; i < argc; i++){
    const std::string a = argv[i];
    auto arg_val = [&](const char* opt_name){
      if(i + 1 >= argc){
        std::cerr << "Error: option '" << opt_name << "' requires an argument.\n";
        std::exit(1);
      }

      i++;
      return argv[i];
    };

    if(a == "-h" || a == "--help"){
      print_help(argv[0]);
      std::exit(0);
    }else if(a == "-m" || a == "--mig"){
      opt.mig_file = arg_val(a.c_str());
    }else if(a == "-d" || a == "--data"){
      opt.baypass_file = arg_val(a.c_str());
    }else if(a == "-s" || a == "--site"){
      int v;
      if (!parse_int_c(arg_val(a.c_str()), v)) {
        std::cerr << "Invalid --site" << std::endl;
        std::exit(1);
      }
      opt.site = v;
    }else if(a == "-o" || a == "--output"){
      opt.simulated_file = arg_val(a.c_str());
    }else if(a == "--pval"){
      opt.pval_file = arg_val(a.c_str());
    }else if(a == "-e" || a == "--env"){
      opt.env_file = arg_val(a.c_str());
    }else if(a == "-l" || a == "--log-file"){
      opt.log_file = arg_val(a.c_str());
    }else if(a == "-b" || a == "--bin-file"){
      opt.bin_file = arg_val(a.c_str());
    }else if(a == "--mode"){
      int v;
      if (!parse_int_c(arg_val(a.c_str()), v) || (v != 0 && v != 1)) {
        std::cerr << "Invalid --mode" << std::endl;
        std::exit(1);
      }
      opt.mode = v;
    }else if(a == "--sim-mode"){
      int v;
      if (!parse_int_c(arg_val(a.c_str()), v) || (v != 0 && v != 1)) {
        std::cerr << "Invalid --sim-mode" << std::endl;
        std::exit(1);
      }
      opt.sim_mode = v;
    }else if(a == "--q-tail"){
      double v;

      if(!parse_double_c(arg_val("--q-tail"), v) || (v <= 0.5 || v >= 1.0)){
        std::cerr << "Invalid --q-tail" << std::endl;
        std::exit(1);
      }
      opt.q_tail = v;
    }else if(a == "-w" || a == "--win-width"){
      double v;

      if(!parse_double_c(arg_val("--win-width"), v) || (v < 0.0 || v > 0.5)){
        std::cerr << "Invalid --win-width (max: 0.5)" << std::endl;
        std::exit(1);
      }
      opt.window = v;
    }else if(a == "--mut-density"){
      double v;

      if(!parse_double_c(arg_val("--mut-density"), v) || (v <= 0.0)){
        std::cerr << "Invalid --mut-density" << std::endl;
        std::exit(1);
      }
      opt.mut_density = v;
    }else if(a == "--prerun"){
      int v;
      if (!parse_int_c(arg_val(a.c_str()), v) || v <= 0) {
        std::cerr << "Invalid --prerun" << std::endl;
        std::exit(1);
      }
      opt.pre_run_tree = v;
    }else if(a == "--seed"){
      unsigned int v;
      if(!parse_uint_c(arg_val("--seed"), v)){
        std::cerr << "Invalid --seed" << std::endl;
        std::exit(1);
      }
      opt.seed = v;
      opt.set_seed = 1;
    }else if(a.rfind("-", 0) == 0){
      std::cerr << "Unknown option: " << a << "\n";
      print_help(argv[0]);
      std::exit(1);
    }else{
      std::cerr << "Error: unexpected positional argument: " << a << "\n";
      print_help(argv[0]);
      std::exit(1);
    }
  }

  if(opt.mig_file.empty()){
    std::cerr << "Error: migration file is required.\n";
    print_help(argv[0]);
    std::exit(1);
  }
  if(opt.baypass_file.empty()){
    std::cerr << "Error: data file is required.\n";
    print_help(argv[0]);
    std::exit(1);
  }
  if(opt.site <= 0){
    std::cerr << "Error: a positive number of sites is required.\n";
    print_help(argv[0]);
    std::exit(1);
  }
  if(opt.mode == 1){
    if(opt.env_file.empty()){
      std::cerr << "Error: env file is required for mode 1.\n";
      print_help(argv[0]);
      std::exit(1);
    }

    if(opt.sim_mode == 1){
      std::cerr << "Error: bin-based analysis is not implemented for GEA mode.\n";
      print_help(argv[0]);
      std::exit(1);
    }
  }
  if(opt.sim_mode == 1){
    if(opt.bin_file.empty()){
      std::cerr << "Error: bin file is required for bin-based simulation mode (sim mode 1).\n";
      print_help(argv[0]);
      std::exit(1);
    }
  }else{
    if(!opt.bin_file.empty()){
      std::cerr << "Warning: bin file is not used for window-based simulation mode (sim mode 0)\n\n";
    }
  }

  opt.version = print_version();

  return opt;
}
