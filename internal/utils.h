/*
 * utils.h
 *
 *  Created on: Mar 31, 2015
 *      Author: nicola
 */

#include <iostream>

using namespace std;

#ifndef UTILS_LZRLCSA_H_
#define UTILS_LZRLCSA_H_

typedef uint32_t uint;
typedef uint64_t ulint;

string get_time(uint64_t time){

	stringstream ss;

	if(time>=3600){

		uint64_t h = time/3600;
		uint64_t m = (time%3600)/60;
		uint64_t s = (time%3600)%60;

		ss  << time << " seconds. ("<< h << "h " << m << "m " << s << "s" << ")";

	}else if (time>=60){

		uint64_t m = time/60;
		uint64_t s = time%60;

		ss << time << " seconds. ("<< m << "m " << s << "s" << ")";

	}else{

		ss << time << " seconds.";

	}

	return ss.str();

}

//parse pizza&chilli patterns header:
void header_error(){
	cout << "Error: malformed header in patterns file" << endl;
	cout << "Take a look here for more info on the file format: http://pizzachili.dcc.uchile.cl/experiments.html" << endl;
	exit(0);
}

ulint get_number_of_patterns(string header){

	ulint start_pos = header.find("number=");
	if (start_pos == std::string::npos or start_pos+7>=header.size())
		header_error();

	start_pos += 7;

	ulint end_pos = header.substr(start_pos).find(" ");
	if (end_pos == std::string::npos)
		header_error();

	ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

	return n;

}

ulint get_patterns_length(string header){

	ulint start_pos = header.find("length=");
	if (start_pos == std::string::npos or start_pos+7>=header.size())
		header_error();

	start_pos += 7;

	ulint end_pos = header.substr(start_pos).find(" ");
	if (end_pos == std::string::npos)
		header_error();

	ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

	return n;

}



#endif /* UTILS_H_ */
