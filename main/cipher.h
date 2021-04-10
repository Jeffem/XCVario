#ifndef CIPHEREV_H
#define CIPHEREV_H
#include <stdlib.h> // abs
#include <string> // std::string
#include <locale> // std::locale, std::toupper
#include "cipher.h"


class Cipher {
public:
	static std::string Encrypt(std::string key, std::string text);
	static std::string Decrypt(std::string key, std::string text);
	static bool checkKeyAHRS();
	static bool init();
private:
	static std::string id();
	static void FormatKey(std::string&, std::string);
	static void FormatEncrypted(std::string&);
	static char ShiftChar(char, int);
	static char ShiftCharBack(char, int);
	static int CalculateDistance(char, char);
};

#endif
