#include "options.hpp"

std::string print_version(){
  std::string ver = "SimGEA migration estimation v1.1.0";
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
  std::cout << "\nUsage: " << prog << " [options] -i <baypass_style_file>\n\n"
    << "Required:\n"
    << "  -i, --input <STR>            BayPass-style SNP count file (space-separated).\n\n"
    << "Optional (defaults in [brackets]):\n"
    << "  -o, --out-prefix <STR>       Output prefix [out]\n"
    << "      --maf <DOUBLE>           MAF filter in [0,0.5) [0.0]\n"
    << "      --lr <DOUBLE>            AdamW learning rate [1e-3]\n"
    << "      --beta1 <DOUBLE>         AdamW beta1 [0.9]\n"
    << "      --beta2 <DOUBLE>         AdamW beta2 [0.999]\n"
    << "      --weight-decay <DOUBLE>  AdamW weight decay [0.01]\n"
    << "      --max-step <INT>         Max optimization steps [100000]\n"
    << "      --min-error <DOUBLE>     Relative improvement threshold [1e-6]\n"
    << "      --min-abs-error <DOUBLE> Absolute improvement threshold [0.01]\n"
    << "      --seed <UNSIGNED INT>    Random seed [randomly determined]\n"
    << "  -l, --log-file <STR>         Log file [log.txt]\n"
    << "  -h, --help                   Show this help and exit\n";
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
      opt.baypass_file = arg_val(a.c_str());
    }else if(a == "-o" || a == "--out-prefix"){
      opt.out_prefix = arg_val(a.c_str());
    }else if(a == "-l" || a == "--log-file"){
      opt.log_file = arg_val(a.c_str());
    }else if(a == "--maf"){
      double v;
      
      if(!parse_double_c(arg_val("--maf"), v)){ 
        std::cerr << "Invalid --maf" << std::endl;
        std::exit(1);
      }
      opt.maf_filter = v;
    }else if(a == "--lr"){
      double v;
      
      if(!parse_double_c(arg_val("--lr"), v)){ 
        std::cerr << "Invalid --lr" << std::endl;
        std::exit(1);
      } 
      opt.lr = v;
    }else if(a == "--beta1"){
      double v;
      
      if(!parse_double_c(arg_val("--beta1"), v)){
        std::cerr << "Invalid --beta1" << std::endl;
        std::exit(1);
      }
      opt.beta1 = v;
    }else if(a == "--beta2"){
      double v;
      if(!parse_double_c(arg_val("--beta2"), v)){
        std::cerr << "Invalid --beta2" << std::endl;
        std::exit(1);
      }
      opt.beta2 = v;
    }else if(a == "--weight-decay"){
      double v; if(!parse_double_c(arg_val("--weight-decay"), v)){
        std::cerr << "Invalid --weight-decay" << std::endl;
        std::exit(1);
      }
      opt.weight_decay = v;
    }else if(a == "--max-step"){
      int v;
      if(!parse_int_c(arg_val("--max-step"), v)){
        std::cerr << "Invalid --max-step" << std::endl;
        std::exit(1);
      }
      opt.max_step = v;
    }else if(a == "--min-error"){
      double v;
      if(!parse_double_c(arg_val("--min-error"), v)){
        std::cerr << "Invalid --min-error" << std::endl;
        std::exit(1);
      }
      opt.min_error = v;
    }else if(a == "--min-abs-error"){
      double v;
      if(!parse_double_c(arg_val("--min-abs-error"), v)){
        std::cerr << "Invalid --min-abs-error" << std::endl;
        std::exit(1);
      }
      opt.min_abs_error = v;
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

  if(opt.baypass_file.empty()){
    std::cerr << "Error: input file is required.\n";
    print_help(argv[0]);
    std::exit(1);
  }

  if(opt.maf_filter < 0.0){
    opt.maf_filter = 0.0;
  }
  if(opt.maf_filter > 0.5){
    std::cerr << "Warning: --maf > 0.5; clamping to 0.5\n";
    opt.maf_filter = 0.5;
  }

  opt.version = print_version();

  return opt;
}
