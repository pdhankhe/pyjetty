#include <iostream>

void testing() {
    std::string str = "hello %s hi";
    std::string full_str = Form(str.c_str(), "you");
    std::cout << full_str << std::endl;  // Output: hello you hi
}