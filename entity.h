
#ifndef _ENTITIES_H_
#define _ENTITIES_H_

#include <string>
#include <vector>

class Entity {
public:
	Entity() {}
	virtual ~Entity() {}

	virtual bool fullEncrypt() const { return true; }

	// For encryption/decryption
	virtual std::string getType() const = 0;
	virtual int getID() const = 0;
	virtual void setID(int new_id) = 0;

	virtual std::vector<int> getValues() const = 0;
	virtual void setValues(const std::vector<int>& values) = 0;

	virtual std::vector<int> getKsValues() const = 0;
	virtual void setKsValues(const std::vector<int>& values) = 0;
};

#endif
