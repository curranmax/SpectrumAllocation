
#include "key_server.h"

#include "entity.h"

#include "cryptopp/aes.h"
#include "cryptopp/modes.h"
#include "cryptopp/osrng.h"

#include <iostream>
#include <map>
#include <stdlib.h>
#include <vector>

using namespace CryptoPP;

void KeyServer::encryptAndInit(Entity* en) {
	if(en->getID() < 0) {
		en->setID(next_id);
		// next_id++;

		// Init key and iv for this entity
		int en_id = en->getID();

		if(aes_vals.count(en_id) == 0) {
			SecByteBlock key(0x00, AES::DEFAULT_KEYLENGTH);
			rnd.GenerateBlock(key, key.size());
	
			SecByteBlock iv(AES::BLOCKSIZE);
			rnd.GenerateBlock(iv, iv.size());
	
			aes_vals[en_id] = new AesPair(key, iv);
		}
	}
	encrypt(en);
}

void KeyServer::encrypt(Entity* en) {
	if(!en->fullEncrypt()) {
		// TODO XOR with some value
		return;
	}

	// pt_vs is the plaintext values of en. en_vs will be the encrypted version of
	std::vector<int> pt_vs = en->getKsValues();
	std::vector<int> en_vs(pt_vs.size());


	std::vector<byte> bytes(4);
	for(unsigned int i = 0; i < pt_vs.size(); ++i) {
		toBytes(pt_vs[i], &bytes);

		auto this_aes_vals = aes_vals.find(en->getID());
		if(this_aes_vals == aes_vals.end()) {
			std::cerr << "No aes vals for: " << en->getID() << std::endl;
			exit(1);
		}

		this_aes_vals->second->encryptor().ProcessData(bytes.data(), bytes.data(), bytes.size());

		int encrypted_v = fromBytes(bytes);
		en_vs[i] = encrypted_v;

		/*{
			// Check that encryption was done right
			std::vector<byte> check_bytes(4);
			toBytes(encrypted_v, &check_bytes);
			
			auto this_aes_vals = aes_vals.find(en->getID());
			if(this_aes_vals == aes_vals.end()) {
				std::cerr << "No aes vals for: " << en->getID() << std::endl;
				exit(1);
			}
			this_aes_vals->second->decryptor().ProcessData(check_bytes.data(), check_bytes.data(), check_bytes.size());

			int check_val = fromBytes(check_bytes);
			
			if(check_val != pt_vs[i]) {
				std::cerr << "Encryption was not symmetric" << std::endl;
				exit(1);
			}
		}*/
	}

	en->setKsValues(en_vs);
}

void KeyServer::decrypt(Entity* en) {
	if(!en->fullEncrypt()) {
		// TODO XOR with some value
		return;
	}

	std::vector<int> en_vs = en->getKsValues();
	std::vector<int> pt_vs(en_vs.size());

	std::vector<byte> bytes(4);
	for(unsigned int i = 0; i < en_vs.size(); ++i) {
		toBytes(en_vs[i], &bytes);

		auto this_aes_vals = aes_vals.find(en->getID());
		if(this_aes_vals == aes_vals.end()) {
			std::cerr << "No aes vals for: " << en->getID() << std::endl;
			exit(1);
		}
		this_aes_vals->second->decryptor().ProcessData(bytes.data(), bytes.data(), bytes.size());

		int decrypted_v = fromBytes(bytes);
		pt_vs[i] = decrypted_v;

		/*{
			// Check that encryption was done right
			std::vector<byte> check_bytes(4);
			toBytes(decrypted_v, &check_bytes);

			auto this_aes_vals = aes_vals.find(en->getID());
			if(this_aes_vals == aes_vals.end()) {
				std::cerr << "No aes vals for: " << en->getID() << std::endl;
				exit(1);
			}
			this_aes_vals->second->encryptor().ProcessData(check_bytes.data(), check_bytes.data(), check_bytes.size());

			int check_val = fromBytes(check_bytes);

			if(check_val != en_vs[i]) {
				std::cerr << "Encryption was not symmetric" << std::endl;
				exit(1);
			}
		}*/
	}

	en->setKsValues(pt_vs);
}

int KeyServer::encryptPRThreshold(const PRint& pr) {
	int pt_thresh = pr.threshold;

	std::vector<byte> bytes(4);
	toBytes(pt_thresh, &bytes);

	auto this_aes_vals = aes_vals.find(pr.getID());
	if(this_aes_vals == aes_vals.end()) {
		std::cerr << "No aes vals for: " << pr.getID() << std::endl;
		exit(1);
	}

	this_aes_vals->second->encryptor().ProcessData(bytes.data(), bytes.data(), bytes.size());

	int en_pt_thresh = fromBytes(bytes);

	return en_pt_thresh;
}

void KeyServer::toBytes(int v, std::vector<byte>* bytes) const {
	for(int i = 0; i < 4; ++i) {
		(*bytes)[i] = (v >> (i * 8));
	}
}

int KeyServer::fromBytes(const std::vector<byte>& bytes) const {
	int v = 0;
	for(int i = 0; i < 4; ++i) {
		v += (int(bytes[i]) << (i * 8));
	}
	return v;
}

