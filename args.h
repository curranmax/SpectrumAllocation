
#ifndef _ARGS_H_
#define _ARGS_H_

#include <string>
#include <vector>

class Args {
  private:
	// Container used to declare command line input
	class Arg_Val {
	  public:
		Arg_Val(float* val_ptr_,float def_val,const std::string &flag_str_,const std::string &input_name_,const std::string &help_text_) { val_ptr = (void*) val_ptr_; *val_ptr_ = def_val; flag_str = flag_str_; val_type = "f"; input_name = input_name_; help_text = help_text_; }
		Arg_Val(int* val_ptr_,int def_val,const std::string &flag_str_,const std::string &input_name_,const std::string &help_text_) { val_ptr = (void*) val_ptr_; *val_ptr_ = def_val; flag_str = flag_str_; val_type = "i"; input_name = input_name_; help_text = help_text_; }
		Arg_Val(std::string* val_ptr_,const std::string &def_val,const std::string &flag_str_,const std::string &input_name_,const std::string &help_text_) { val_ptr = (void*) val_ptr_; *val_ptr_ = def_val; flag_str = flag_str_; val_type = "s"; input_name = input_name_; help_text = help_text_; }
		Arg_Val(bool* val_ptr_,const std::string &flag_str_,const std::string &help_text_) { val_ptr = (void*) val_ptr_; *val_ptr_ = false; flag_str = flag_str_; val_type = "b"; help_text = help_text_; }

		Arg_Val(const Arg_Val& av) {val_ptr = av.val_ptr; flag_str = av.flag_str; val_type = av.val_type; input_name = av.input_name; help_text = av.help_text; }

		bool isMatch(const std::vector<std::string> &cmd_args, unsigned int &i);
		void setVal(const std::vector<std::string> &cmd_args, unsigned int &i);

		void* val_ptr;
		std::string flag_str;
		std::string val_type;
		std::string input_name;
		std::string help_text;
	};
	
  public:
	Args(const int &argc, char const *argv[]);
	~Args() {}

	void printHelp() const;
	std::vector<Arg_Val> args;

	// General
	int rand_seed;
	bool brief_out;

	// Generation params
	int num_pu;
	int num_ss;
	int num_su;
	float location_range;

	// Propagation Model
	std::string propagation_model;

	// Log Distance Params
	float ld_path_loss0; // in dBm
	float ld_dist0;
	float ld_gamma;

	// Spectrum Manager params
	int num_ss_selection;

	// TODO allow this to be a floating point value
	int ss_receive_power_alpha;
	int path_loss_alpha;

	// The number of bits to shift ints in the S2-PC calculations.
	int num_float_bits;

	// Number of bits overall to use in S2-PC.
	int s2_pc_bit_count;
};

#endif
