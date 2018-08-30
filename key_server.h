
#ifndef _KEY_SERVER_H_
#define _KEY_SERVER_H_

#include "entity.h"
#include "primary_user.h"

#include "cryptopp/aes.h"
#include "cryptopp/modes.h"
#include "cryptopp/osrng.h"

#include <map>
#include <vector>

class AesPair {
public:
	AesPair(CryptoPP::SecByteBlock _key, CryptoPP::SecByteBlock _iv) : key(_key), iv(_iv) {}

	AesPair() = delete;
	~AesPair() {}
	AesPair(const AesPair& aes_pair) = delete;

	const AesPair& operator=(const AesPair& aes_pair) = delete;

	CryptoPP::CFB_Mode<CryptoPP::AES>::Encryption encryptor() const {
		return CryptoPP::CFB_Mode<CryptoPP::AES>::Encryption(key, key.size(), iv);
	}

	CryptoPP::CFB_Mode<CryptoPP::AES>::Decryption decryptor() const {
		return CryptoPP::CFB_Mode<CryptoPP::AES>::Decryption(key, key.size(), iv);
	}

	CryptoPP::SecByteBlock key, iv;
};

class KeyServer {
public:
	KeyServer() : rnd(), next_id(0), aes_vals() {}
	~KeyServer() {
		for(auto itr = aes_vals.begin(); itr != aes_vals.end(); ++itr) {
			delete itr->second;
		}
	}

	KeyServer(const KeyServer& ks) = delete;
	const KeyServer& operator=(const KeyServer& ks) = delete;

	void encryptAndInit(Entity* en);
	void encrypt(Entity* en);
	void decrypt(Entity* en);

	int encryptPRThreshold(const PRint& pr);
private:
	void toBytes(int v, std::vector<CryptoPP::byte>* bytes) const;
	int fromBytes(const std::vector<CryptoPP::byte>& bytes) const;

	CryptoPP::AutoSeededRandomPool rnd;

	int next_id;
	std::map<int, AesPair*> aes_vals;
};

#endif
