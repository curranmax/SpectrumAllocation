
#ifndef _SU_SERVER_H_
#define _SU_SERVER_H_

#include "secondary_user.h"

#include <vector>

class SUServer {
public:
	SUServer() = delete;
	SUServer(const std::vector<SU>& _sus) : sus(_sus) {}

	SUServer(const SUServer& su_server) = delete;
	const SUServer& operator=(const SUServer& su_server) = delete;

	void run() const;

private:
	bool anyUnused(const std::vector<bool>& used) const;

	std::vector<SU> sus;
};

#endif
