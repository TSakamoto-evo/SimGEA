#include "options.hpp"

std::string print_version(){
  std::string ver = "under development v0.0";
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


void print_help(const char* prog){
  print_version();
  std::cout << "\nUsage: " << prog << " [options] -i <migration_file> -s <site_number> -a <allele_number> \n\n"

    << "Required:\n"
    << "  -i, --input <STR>           Migration file (space-separated).\n"
    << "  -s  --site <INT>            Number of sites in the output file.\n"
    << "  -a  --allele <INT>          Number of alleles taken from each subpopulation.\n\n"
    << "Optional (defaults in [brackets]):\n"
    << "  -o, --output <STR>          Output baypass-styled file [out_data.txt]\n"
    << "      --n-pre <INT>           Number of trees in the prerun [10000]\n"
    << "      --maf <DOUBLE>          MAF filter in [0,0.5) [0.0]\n"
    << "      --mut-density <DOUBLE>  Expected number of mutations per tree [0.25]\n"
    << "  -h, --help                  Show this help and exit\n";
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
    }else if(a == "-i" || a == "--input"){
      opt.mig_file = arg_val(a.c_str());
    }else if(a == "-o" || a == "--output"){
      opt.baypass_file = arg_val(a.c_str());
    }else if(a == "-s" || a == "--site"){
      int v;
      if (!parse_int_c(arg_val(a.c_str()), v)) {
        std::cerr << "Invalid --site" << std::endl;
        std::exit(1);
      }
      opt.site = v;
    }else if(a == "-a" || a == "--allele"){
      int v;
      if (!parse_int_c(arg_val(a.c_str()), v)) {
        std::cerr << "Invalid --allele" << std::endl;
        std::exit(1);
      }
      opt.sample_num = v;
    }else if(a == "--maf"){
      double v;

      if(!parse_double_c(arg_val("--maf"), v)){
        std::cerr << "Invalid --maf" << std::endl;
        std::exit(1);
      }
      opt.freq_cut = v;
    }else if(a == "--n-pre"){
      int v;
      if (!parse_int_c(arg_val(a.c_str()), v)) {
        std::cerr << "Invalid --n-pre" << std::endl;
        std::exit(1);
      }
      opt.pre_run_tree = v;
    }else if(a == "--mut-density"){
      double v;

      if(!parse_double_c(arg_val("--mut-density"), v)){
        std::cerr << "Invalid --mut-density" << std::endl;
        std::exit(1);
      }
      opt.exp_mut_per_tree = v;
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
    std::cerr << "Error: input file is required.\n";
    print_help(argv[0]);
    std::exit(1);
  }

  if(opt.site <= 0){
    std::cerr << "Error: positive site number is required.\n";
    print_help(argv[0]);
    std::exit(1);
  }
  if(opt.sample_num <= 0){
    std::cerr << "Error: positive sample number is required.\n";
    print_help(argv[0]);
    std::exit(1);
  }

  if(opt.freq_cut < 0.0){
    opt.freq_cut = 0.0;
  }
  if(opt.freq_cut > 0.5){
    std::cerr << "Warning: --maf > 0.5; clamping to 0.5\n";
    opt.freq_cut = 0.5;
  }

  opt.version = print_version();

  return opt;
}
