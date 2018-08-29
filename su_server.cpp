
#include "su_server.h"

#include "split.h"

#include "cryptoTools/Network/IOService.h"

#include <iostream>

using namespace osuCrypto;

void SUServer::run() const {
	// Connect to both SMs
	IOService ios(1);
	Endpoint ep0(ios, "127.0.0.1:1213", EpMode::Server, "sm0-su");
	Channel ch0 = ep0.addChannel("sm0-su-channel");
	ch0.waitForConnection();

	Endpoint ep1(ios, "127.0.0.1:1214", EpMode::Server, "sm1-su");
	Channel ch1 = ep1.addChannel("sm1-su-channel");
	ch1.waitForConnection();

	// Wait for message from the SMs
	std::vector<bool> used_su(sus.size(), false);

	unsigned int x = 0;
	while(anyUnused(used_su)) {
		if(x >= sus.size()) {
			std::cerr << "Expected to be finished by now" << std::endl;
			exit(1);
		}

		std::array<int, 3> sm0_vals;
		std::array<int, 3> sm1_vals;

		ch0.recv(sm0_vals);
		ch1.recv(sm1_vals);

		if(sm0_vals[0] != sm1_vals[0]) {
			std::cerr << "Mismatched indexes: " << sm0_vals[0] << ", " << sm1_vals[0] << std::endl;
			exit(1);
		}
		unsigned int this_index = sm0_vals[0];

		if(this_index < 0 || this_index >= sus.size()) {
			std::cerr << "Invalid index: " << this_index << std::endl;
			exit(1);
		}

		if(used_su[this_index]) {
			std::cerr << "Repeat index: " << this_index << std::endl;
		}

		used_su[this_index] = true;

		if(sm0_vals[2] != sm1_vals[2]) {
			std::cerr << "Mismatched factor: " << sm0_vals[2] << ", " << sm1_vals[2] << std::endl;
			exit(1);
		}
		float factor = sm0_vals[2];

		float max_tp = float(sm0_vals[1] + sm1_vals[1]) / factor;

		float this_tp = max_tp + sus[this_index].less_max_tp;
		int this_transmit = (this_tp > sus[this_index].min_tp ? 1 : 0);

		auto this_tp_split = splitInt(int(this_tp * factor));
		auto this_transmit_split = splitInt(this_transmit);

		ch0.send(std::array<int, 2>{this_tp_split.first,  this_transmit_split.first});
		ch1.send(std::array<int, 2>{this_tp_split.second, this_transmit_split.second});

		++x;
	}

	// Compute the transmit power and if it will transmit

}

bool SUServer::anyUnused(const std::vector<bool>& used) const {
	for(unsigned int i = 0; i < used.size(); ++i) {
		if(!used[i]) {
			return true;
		}
	}
	return false;
}
