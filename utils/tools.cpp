/***************************************************************************
 *   Copyright (C) 2009-2015 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
 *                                                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/



#include "tools.h"
#include "starttree.h" //for START_TREE_RECOGNIZED macro.
#include "timeutil.h"
#include "MPIHelper.h"
#ifndef CLANG_UNDER_VS
    #include <dirent.h>
#else
     //James B. Workaround for Windows builds where these macros might not be defined
    #ifndef S_ISDIR
    #define S_ISDIR(mode) (((mode) & S_IFMT) == S_IFDIR)
    #endif
    #ifndef S_ISREG
    #define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
    #endif
#endif
#include <thread>


#if defined(Backtrace_FOUND)
#include <execinfo.h>
#include <cxxabi.h>
#endif

#include "tools.h"
#include "timeutil.h"
#include "progress.h"
#include "gzstream.h"
#include "MPIHelper.h"
#include "alignment/alignment.h"

VerboseMode verbose_mode;
extern void printCopyright(ostream &out);

#if defined(WIN32)
#include <sstream>
#endif

/********************************************************
        Miscellaneous
 ********************************************************/

/**
        Output an error to screen, then exit program
        @param error error message
 */
void outError(const char *error, bool quit) {
	if (error == ERR_NO_MEMORY) {
        print_stacktrace(cerr);
	}
#ifndef BUILD_LIB
	cerr << error << endl;
#endif
    if (quit) {
#ifndef BUILD_LIB
        exit(2);
#else
        throw runtime_error(error);
#endif
    }
}

/**
        Output an error to screen, then exit program
        @param error error message
 */
void outError(string error, bool quit) {
    outError(error.c_str(), quit);
}

void outError(const char *error, const char *msg, bool quit) {
    string str = error;
    str += msg;
    outError(str, quit);
}

void outError(const char *error, string msg, bool quit) {
    string str = error;
    str += msg;
    outError(str, quit);
}

/**
        Output a warning message to screen
        @param error warning message
 */
void outWarning(const char *warn) {
    cout << "WARNING: " << warn << endl;
}

void outWarning(string warn) {
    outWarning(warn.c_str());
}

double tryGeneratingBlength(Params &params) {
    // randomly generate branch lengths based on
    // a user-specified distribution
    if (params.branch_distribution)
        return random_number_from_distribution(params.branch_distribution, true);
    // or an exponential distribution (by default)
    else
    {
        double len = random_double_exponential_distribution(params.mean_len);

        if (len < params.min_len) {
            int fac = random_int(1000);
            double delta = static_cast<double> (fac) / 1000.0; //delta < 1.0
            len = params.min_len + delta / 1000.0;
        }

        if (len > params.max_len) {
            int fac = random_int(1000);
            double delta = static_cast<double> (fac) / 1000.0; //delta < 1.0
            len = params.max_len - delta / 1000.0;
        }
        return len;
    }
}

double randomLen(Params &params) {
    // bug fixed: avoid negative branch lengths
    double len = -1;
    int attemp = 0;
    while ((len < params.min_len || len > params.max_len) && attemp < 1000)
    {
        len = tryGeneratingBlength(params);
        ++attemp;
    }
    if (len < params.min_len || len > params.max_len)
        outError("Failed to generate a branch length (in the range(" + convertDoubleToString(params.min_len) + ", " + convertDoubleToString(params.max_len) + ")) after 1000 attempts. Please check the input and try again!");
    return len;
}

std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

//From Tung

string convertIntToString(int number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

string convertInt64ToString(int64_t number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

string convertDoubleToString(double number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

bool iEquals(const string a, const string b)
{
    unsigned int sz = a.size();
    if (b.size() != sz)
        return false;
    for (unsigned int i = 0; i < sz; ++i)
        if (tolower(a[i]) != tolower(b[i]))
            return false;
    return true;
}

//From Tung

bool copyFile(const char SRC[], const char DEST[]) {
    std::ifstream src; // the source file
    std::ofstream dest; // the destination file

    src.open(SRC, std::ios::binary); // open in binary to prevent jargon at the end of the buffer
    dest.open(DEST, std::ios::binary); // same again, binary
    if (!src.is_open() || !dest.is_open())
        return false; // could not be copied

    dest << src.rdbuf(); // copy the content
    dest.close(); // close destination file
    src.close(); // close source file

    return true; // file copied successfully
}

bool fileExists(string strFilename) {
    struct stat stFileInfo;
    bool blnReturn;
    int intStat;

    // Attempt to get the file attributes
    intStat = stat(strFilename.c_str(), &stFileInfo);
    if (intStat == 0) {
        // We were able to get the file attributes
        // so the file obviously exists.
        blnReturn = true;
    } else {
        // We were not able to get the file attributes.
        // This may mean that we don't have permission to
        // access the folder which contains this file. If you
        // need to do that level of checking, lookup the
        // return values of stat which will give you
        // more details on why stat failed.
        blnReturn = false;
    }
    return (blnReturn);
}

int isDirectory(const char *path) {
    struct stat statbuf;
    if (stat(path, &statbuf) != 0)
        return 0;
    return S_ISDIR(statbuf.st_mode);
}

int isFile(const char *path) {
    struct stat statbuf;
    if (stat(path, &statbuf) != 0)
        return 0;
    return S_ISREG(statbuf.st_mode);
}

int getFilesInDir(const char *path, StrVector &filenames)
{
    if (!isDirectory(path))
        return 0;
    string path_name = path;
    if (path_name.back() != '/')
        path_name.append("/");
#ifndef CLANG_UNDER_VS
    DIR *dp;
    struct dirent *ep;
    dp = opendir (path);
    
    if (dp != nullptr)
    {
        while ((ep = readdir (dp)) != nullptr) {
            if (isFile((path_name + ep->d_name).c_str()))
                filenames.push_back(ep->d_name);
        }
        
        (void) closedir (dp);
        return 1;
    }
    else
        return 0;
    
    return 1;
#else
    //TODO: FindFirstFile et cetera, et cetera, et cetera
    return 0;
#endif
}
int convert_int(const char *str) {
    char *endptr;
    int i = strtol(str, &endptr, 10);

    if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        outError(err);
    }

    return i;
}

int convert_int(const char *str, int &end_pos) {
	char *endptr;
	int i = strtol(str, &endptr, 10);

	if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL) {
		string err = "Expecting integer, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}
	end_pos = endptr - str;
	return i;
}

void convert_int_vec(const char *str, IntVector &vec) {
    char *beginptr = (char*)str, *endptr;
    vec.clear();
    do {
		int i = strtol(beginptr, &endptr, 10);

		if ((i == 0 && endptr == beginptr) || abs(i) == HUGE_VALL) {
			string err = "Expecting integer, but found \"";
			err += beginptr;
			err += "\" instead";
			throw err;
		}
		vec.push_back(i);
		if (*endptr == ',') endptr++;
		beginptr = endptr;
    } while (*endptr != 0);
}


int64_t convert_int64(const char *str) {
    char *endptr;
    int64_t i = (int64_t)strtoll(str, &endptr, 10); // casted because 'long long' may be larger than int64_t

    if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting large integer , but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }

    return i;
}

int64_t convert_int64(const char *str, int &end_pos) {
	char *endptr;
	int64_t i = (int64_t)strtoll(str, &endptr, 10); // casted because 'long long' may be larger than int64_t

	if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL) {
		string err = "Expecting large integer, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}
	end_pos = endptr - str;
	return i;
}


double convert_double(const char *str) {
    char *endptr;
    double d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || *endptr != 0) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    return d;
}

double convert_double(const char *str, int &end_pos) {
	char *endptr;
	double d = strtod(str, &endptr);
	if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF) {
		string err = "Expecting floating-point number, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}
	end_pos = endptr - str;
	return d;
}

double convert_double_with_distribution(const char *str, int &end_pos, bool non_negative, char separator)
{
    // convert normal double
    char *endptr;
    double d = strtod(str, &endptr);
    
    // generate a double from a distribution
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF) {
        string tmp_str(str);
        size_t pos = tmp_str.find(separator);
        if (pos!= std::string::npos)
            end_pos = pos;
        else
            end_pos = tmp_str.length();
            
        d = random_number_from_distribution(tmp_str.substr(0, end_pos), non_negative);
    }
    else
        end_pos = endptr - str;
    return d;
}

void convert_double_vec(const char *str, DoubleVector &vec, char separator) {
    char *beginptr = (char*)str, *endptr;
    vec.clear();
    do {
		double d = strtod(beginptr, &endptr);

		if ((d == 0.0 && endptr == beginptr) || fabs(d) == HUGE_VALF) {
			string err = "Expecting floating-point number, but found \"";
			err += beginptr;
			err += "\" instead";
			throw err;
		}
		vec.push_back(d);
		if (*endptr == separator) endptr++;
		beginptr = endptr;
    } while (*endptr != 0);
}

void convert_double_vec_with_distributions(const char *str, DoubleVector &vec, bool non_negative, char separator)
{
    string tmp_str(str);
    vec.clear();
    
    // extract/generate double numbers one by one
    while (tmp_str.length() > 0) {
        // extract sub-string by separator
        size_t pos = tmp_str.find(separator);
        string token = tmp_str.substr(0, pos);
        
        // convert/generate a double
        double d = convert_double_with_distribution(token.c_str(), non_negative);
        vec.push_back(d);
        
        // remove the current double/distribution name from tmp_str
        if(pos != std::string::npos)
            tmp_str.erase(0, pos + 1);
        else
            tmp_str = "";
    }
}

void convert_double_array_with_distributions(string tmp_str, double* array, int num_items, bool non_negative, char separator)
{
    // validate the number of items
    size_t num_separators = std::count(tmp_str.begin(), tmp_str.end(), separator);
    if (num_separators != num_items - 1)
        outError("The number of items in "+tmp_str+" is "+convertIntToString(num_separators+1)+", which is different from the expected number of items ("+convertIntToString(num_items)+"). Please check and try again!");

    // extract/generate double numbers one by one
    for (int i = 0; i < num_items; i++) {
        // extract sub-string by separator
        size_t pos = tmp_str.find(separator);
        string token = tmp_str.substr(0, pos);
        
        // convert/generate a double
        array[i] = convert_double_with_distribution(token.c_str(), non_negative);
        
        // remove the current double/distribution name from tmp_str
        tmp_str.erase(0, pos + 1);
    }
}

string convert_time(const double sec) {
    int sec_int = (int) floor(sec);
    int secs = sec_int % 60;
    int mins = (sec_int % 3600) / 60;
    int hours = sec_int / 3600;
    stringstream ss;
    ss << hours << "h:" << mins << "m:" << secs << "s";
    return ss.str();
}

void convert_range(const char *str, int &lower, int &upper, int &step_size) {
    char *endptr;

    // parse the lower bound of the range
    int d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    //lower = d;
    int d_save = d;
    upper = d;
    if (*endptr == 0) return;


    // parse the upper bound of the range
    str = endptr + 1;
    d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }

    lower = d_save;
    upper = d;
    if (*endptr == 0) return;

    // parse the step size of the range
    str = endptr + 1;
    d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    step_size = d;
}

void convert_range(const char *str, double &lower, double &upper, double &step_size) {
    char *endptr;

    // parse the lower bound of the range
    double d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    //lower = d;
    double d_save = d;
    upper = d;
    if (*endptr == 0) return;


    // parse the upper bound of the range
    str = endptr + 1;
    d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }

    lower = d_save;
    upper = d;
    if (*endptr == 0) return;

    // parse the step size of the range
    str = endptr + 1;
    d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || *endptr != 0) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    step_size = d;
}

void convert_string_vec(const char *str, StrVector &vec, char separator) {
    char *beginptr = (char*)str, *endptr;
    vec.clear();
    string elem;
    do {
    	endptr = strchr(beginptr, separator);
    	if (!endptr) {
    		elem.assign(beginptr);
    		vec.push_back(elem);
    		return;
    	}
    	elem.assign(beginptr, endptr-beginptr);
    	vec.push_back(elem);
		beginptr = endptr+1;
    } while (*endptr != 0);

}

bool renameString(string &name) {
    bool renamed = false;
    for (string::iterator i = name.begin(); i != name.end(); i++) {
        if ((*i) == '/') {
            // PLL does not accept '/' in names, turn it off
            if (Params::getInstance().start_tree == STT_PLL_PARSIMONY)
                Params::getInstance().start_tree = STT_PARSIMONY;
        }
        if (!isalnum(*i) && (*i) != '_' && (*i) != '-' && (*i) != '.' && (*i) != '|' && (*i) != '/') {
            (*i) = '_';
            renamed = true;
        }
    }
    return renamed;
}

void read_distributions(char* filepath)
{
    const char* builtin_distributions = R"(Generalized_logistic 0.363799 0.313203 0.277533 0.24235 0.260252 0.321134 0.299891 0.315519 0.269172 0.258165 0.287641 0.309442 0.264017 0.23103 0.178778 0.200087 0.336534 0.337547 0.325379 0.335034 0.306336 0.359459 0.249315 0.388073 0.296979 0.345694 0.338733 0.305404 0.294181 0.303477 0.257679 0.417313 0.290922 0.301826 0.292324 0.179902 0.122071 0.348381 0.33887 0.228999 0.377297 0.296036 0.044523 0.262098 0.295087 0.311826 0.309402 0.302522 0.244771 0.325932 0.256499 0.220488 0.230443 0.241669 0.345664 0.262144 0.267938 0.293822 0.326857 0.376613 0.24789 0.303307 0.332711 0.291757 0.234466 0.380376 0.268955 0.278675 0.290531 0.23741 0.066827 0.252364 0.277417 0.325834 0.299602 0.24487 0.332175 0.300581 0.291123 0.312748 0.078066 0.298195 0.380702 0.190509 0.358828 0.379174 0.337185 0.167293 0.271394 0.201814 0.222829 0.332488 0.215259 0.298748 0.257805 0.290615 0.102878 0.292949 0.298971 0.336142 0.223446 0.329779 0.236269 0.280042 0.34049 0.302792 0.300794 0.234272 0.275786 0.231749 0.257803 0.248692 0.286229 0.305252 0.25298 0.274867 0.278125 0.318123 0.222347 0.260319 0.310544 0.262839 0.289586 0.296508 0.311839 0.261722 0.327792 0.314063 0.308725 0.269523 0.245446 0.298576 0.311922 0.289018 0.361806 0.36013 0.242279 0.349005 0.28352 0.277942 0.27051 0.282086 0.292739 0.261955 0.31282 0.194691 0.325831 0.352321 0.330043 0.30416 0.301521 0.324624 0.198778 0.320556 0.440894 0.332368 0.227977 0.322864 0.282516 0.285885 0.231099 0.290483 0.255771 0.306594 0.313493 0.22475 0.21297 0.27572 0.210148 0.30656 0.367425 0.302816 0.255029 0.317457 0.148676 0.234732 0.309629 0.353117 0.414357 0.317483 0.244311 0.310266 0.327092 0.322577 0.330368 0.254572 0.323599 0.266023 0.330568 0.278868 0.352982 0.31807 0.254994 0.366833 0.252355 0.336477 0.363041 0.34216 0.227844 0.340953 0.298944 0.323224 0.275477 0.247568 0.261933 0.306794 0.26298 0.30456 0.256401 0.327072 0.296593 0.277926 0.306128 0.305057 0.226235 0.261324 0.386823 0.312633 0.334579 0.300794 0.226765 0.258937 0.21015 0.313397 0.318114 0.240953 0.325361 0.23713 0.236193 0.364707 0.09543 0.263138 0.255472 0.209969 0.19782 0.221003 0.306613 0.298886 0.267442 0.286558 0.130796 0.280113 0.282211 0.355578 0.25584 0.289949 0.268929 0.324786 0.250406 0.269572 0.303883 0.268379 0.171748 0.30994 0.345134 0.254873 0.129833 0.34198 0.253411 0.284532 0.278311 0.293196 0.305821 0.281977 0.301114 0.266777 0.304058 0.354363 0.29998 0.312962 0.242894 0.358988 0.309335 0.26064 0.177357 0.280343 0.336571 0.384738 0.261577 0.229671 0.194292 0.284024 0.349357 0.406105 0.277227 0.318222 0.328541 0.234773 0.256208 0.298503 0.267952 0.329015 0.332247 0.256429 0.322977 0.303206 0.041074 0.094734 0.234724 0.262491 0.319019 0.272767 0.325977 0.318982 0.310528 0.322474 0.377355 0.368872 0.212661 0.294858 0.11039 0.337951 0.200809 0.245518 0.281807 0.245854 0.277848 0.376086 0.275291 0.203582 0.318731 0.194965 0.255722 0.319838 0.341335 0.329492 0.235443 0.288418 0.310715 0.301169 0.27953 0.31977 0.201715 0.3581 0.332325 0.190273 0.383453 0.294277 0.313312 0.317989 0.348503 0.309744 0.179324 0.344893 0.340539 0.271633 0.223929 0.122569 0.316808 0.258686 0.367824 0.350559 0.184925 0.287621 0.291921 0.206303 0.275647 0.317818 0.190351 0.317876 0.333233 0.260405 0.314324 0.272614 0.274124 0.229243 0.270435 0.253669 0.282295 0.240033 0.273747 0.198832 0.32967 0.337402 0.209729 0.255561 0.374547 0.292157 0.280284 0.253114 0.326909 0.341647 0.283922 0.328876 0.260656 0.284541 0.326051 0.346184 0.217852 0.311937 0.231578 0.325095 0.326498 0.336968 0.297321 0.281219 0.127107 0.324204 0.352182 0.298148 0.256763 0.266463 0.28917 0.326688 0.268065 0.264716 0.274123 0.356898 0.313407 0.304807 0.327584 0.317122 0.29966 0.054486 0.306675 0.37056 0.302364 0.295535 0.232633 0.351527 0.276551 0.283668 0.234175 0.168121 0.221882 0.217143 0.262378 0.162796 0.286549 0.322606 0.327978 0.166185 0.107094 0.240343 0.331932 0.287218 0.297551 0.226403 0.410376 0.292106 0.306494 0.246621 0.252907 0.279642 0.289623 0.352538 0.337796 0.25536 0.256646 0.309608 0.282115 0.3033 0.307261 0.275992 0.282713 0.2838 0.353441 0.226078 0.292404 0.282489 0.358105 0.281517 0.316241 0.253381 0.292601 0.302286 0.193929 0.251386 0.314674 0.30573 -0.206354 0.184509 0.256875 0.265588 0.285296 0.276962 0.259238 0.216939 0.265915 0.315037 0.274147 0.300123 0.402781 0.262571 0.222779 0.348932 0.219679 0.230462 0.272442 0.230498 0.336596 0.275269 0.31521 0.30337 0.266202 0.198975 0.297142 0.350524 0.291291 0.290985 0.34113 0.22364 0.241448 0.269696 0.328101 0.279475 0.24356 0.376772 0.300235 0.324975 0.259819 0.256142 0.290844 0.240551 0.251107 0.279825 0.303731 0.285856 0.257266 0.309779 0.268429 0.253835 0.327948 0.340819 0.209707 0.346565 0.330817 0.273283 0.316428 0.34028 0.27082 0.19219 0.200584 0.324345 0.301378 0.250436 0.276101 0.352346 0.305654 0.208992 0.284911 0.292299 0.284243 0.431182 0.277142 0.262104 0.34329 0.31014 0.206567 0.221605 0.351522 0.309254 0.328483 0.26427 0.349916 0.287899 0.35866 0.318789 0.24227 0.273739 0.313341 0.300235 0.275819 0.326983 0.268149 0.338917 0.2958 0.356719 0.369453 0.290371 0.303753 0.084036 0.363659 0.277632 0.40944 0.091232 0.371642 0.307421 0.257834 0.328097 0.269456 0.234288 0.288251 0.290446 0.326906 0.235733 0.338238 0.232989 0.337715 0.319809 0.306349 0.346429 0.311162 0.25513 0.300376 0.213724 0.223544 0.342333 0.233747 0.23793 0.195557 0.234264 0.343337 0.327066 0.243337 0.283982 0.25023 0.383851 0.316573 0.260721 0.2365 0.336424 0.25442 0.363219 0.422875 0.292284 0.465584 0.260122 0.310313 0.280185 0.233588 0.244522 0.256319 0.290408 0.2534 0.260705 0.344413 0.382906 0.230257 0.305927 0.29825 0.19388 0.217396 0.268437 0.317264 0.323392 0.258615 0.166411 0.231721 0.113373 0.294681 0.251311 0.327711 0.22927 0.194339 0.239293 0.334808 0.270195 0.298885 0.26811 0.351372 0.370811 0.306145 0.190562 0.280724 0.33142 0.27318 0.232841 0.130027 0.323117 0.326668 0.303674 0.226465 0.173001 0.360679 0.26009 0.328231 0.283678 0.200057 0.278999 0.354285 0.374157 0.371166 0.303174 0.334112 0.306962 0.332735 0.299995 0.313238 0.242983 0.331118 0.314005 0.225765 0.338705 0.272283 0.312879 0.299343 0.205338 0.244927 0.318327 0.283998 0.307697 0.389957 0.327577 0.217277 0.295726 0.146595 0.306907 0.281563 0.264 0.352806 0.120943 0.243214 0.308845 0.255728 0.197542 0.290011 0.353052 0.334456 0.289128 0.201144 0.220763 0.382263 0.264179 0.236047 0.368475 0.213221 0.283559 0.235322 0.261908 0.200954 0.011635 0.251288 0.267163 0.326864 0.27198 0.323665 0.299936 0.125222 0.275398 0.312234 0.268429 0.200154 0.185819 0.293116 0.299084 0.258933 0.283395 0.248175 0.183176 0.344321 0.297136 0.344131 0.173493 0.244146 0.247723 0.257914 0.139585 0.226581 0.188051 0.26786 0.304243 0.2612 0.269216 0.293813 0.33594 0.340036 0.294782 0.136761 0.338027 0.334017 0.284144 0.377558 0.253855 0.239247 0.341585 0.301219 0.222854 0.204955 0.385188 0.303908 0.368828 0.204491 0.373765 0.26339 0.279725 0.258224 0.307475 0.276247 0.335483 0.13706 0.286348 0.314757 0.300761 0.231982 0.29167 0.241961 0.345141 0.289442 0.340516 0.282426 0.270503 0.283341 0.33909 0.320187 0.323418 0.318151 0.255541 0.228704 0.335369 0.297904 0.358993 0.344302 0.206315 0.273064 0.308217 0.125667 0.330191 0.386132 0.244322 0.328132 0.348462 0.28154 0.285001 0.285568 0.212048 0.263071 0.183724 0.145403 0.119225 0.304215 0.298517 0.192418 0.285471 0.257938 0.333621 0.291339 0.244867 0.260013 0.313095 0.305021 0.341727 0.269794 0.144422 0.084487 0.247173 0.29659 0.306122 0.179697 0.354462 0.279102 0.357901 0.297665 0.31799 0.244045 0.350837 0.290353 0.292364 0.344325 0.251891 0.308069 0.248317 0.316096 0.222557 0.259914 0.370477 0.179602 0.336114 0.304706 0.266185 0.317406 0.321365 0.336165 0.32789 0.279106 0.320172 0.355546 0.242814 0.13552 0.30337 0.250813 0.237047 0.320047 0.323985 0.347776 0.402441 0.278457 0.307687 0.26086 0.417576 0.319724 0.249277 0.295823 0.360911 0.308001 0.294723 0.277671 0.29898 0.218064 0.251803 0.249479 0.243514 0.24833 0.249269 0.277469 0.366125 0.301805 0.331017 0.259064 0.328109 0.365054 0.206247 0.245856 0.172884 0.313686 0.247422 0.243916 0.35395 0.220123 0.313185 0.254103 0.182547 0.298309 0.289123 0.251778 0.251973 0.318414 0.270089 0.235271 0.249376 0.33918 0.351909 0.327058 0.313975 0.277684 0.330549 0.249332 0.301396 0.337183 0.103806 0.212153 0.359464 0.333055 0.239178 0.336958 0.399539 0.245594 0.339767 0.330008 0.200138 0.272999 0.282139 0.279513 0.334279 0.235428 0.207892 0.21283 0.309937 0.247368 0.187385 0.201026 0.25524 0.353723 0.324293 0.353879 0.431379 0.337694 0.34454 0.355 0.352871 0.302964 0.206905 0.244987 0.312907 0.287111 0.343188 0.269969 0.314492 0.32263 0.331015 0.235652 0.259456 0.300321 0.436108 0.347747 0.325582 0.326411 0.263365 0.336122 0.340153 0.014711 0.28613 0.389232 0.297337 0.331345 0.314414 0.340354 0.320492 0.31058 0.341076 0.324648 0.363326 0.30805 0.223235 0.279956 0.252912 0.312391 0.285699 0.274296 0.31785 0.282243 0.303535 0.3685 0.33568 0.264435 0.27515 0.299322 0.249723 0.35219 0.250568 0.265345 0.243514 0.295327 0.292164 0.317382 0.218889 0.335726 0.301339 0.255064 0.229266 0.433796 0.308085 0.253114 0.277187 0.261273 0.349387 0.359527 0.176964 0.315956 0.283469 0.278288 0.321365 0.212062 0.356948 0.201731 0.216876 0.305531 0.250832 0.365414 0.361186 0.283809 0.384347 0.280245 0.34793 0.406597 0.258394 0.325739 0.36723 0.263813 0.322964 0.319139 0.225419 0.329109 0.355581 0.140803 0.285786 0.263977 0.292095 0.247361 0.249835 0.278165 0.288542 0.260634 0.217882 0.188384 0.311364 0.135016 0.225019 0.328294 0.285733 0.237516 0.198389 0.313495 0.275812 0.19774 0.306504 0.334237 0.226039 0.377486 0.326434 0.431388 0.309237 0.25677 0.178301 0.326125 0.31591 0.247916 0.326454 0.237216 0.276075 0.27966 0.27385 0.254524 0.324365 0.385129 0.296657 0.349485 0.336665 0.378804 0.329152 0.267426 0.259268 0.317341 0.268414 0.272943 0.171274 0.294049 0.219838 0.218484 0.299107 0.216462 0.290662 0.243955 0.265274 0.317071 0.246434 0.243026 0.275551 0.207322 0.193639 0.294482 0.332651 0.339402 0.277366 0.300209 0.263703 0.334661 0.285752 0.260173 0.365594 0.285002 0.389468 0.352159 0.244804 0.350031 0.331839 0.356268 0.290662 0.313201 0.177531 0.32309 0.33779 0.177402 0.253809 0.281584 0.318826 0.332877 0.18783 0.311989 0.351954 0.301753 0.12937 0.36282 0.257102 0.326068 0.315294 0.227995 0.326653 0.206674 0.271749 0.354894 0.265171 0.326021 0.238089 0.329939 0.295873 0.265234 0.236883 0.303479 0.328569 0.216826 0.287449 0.303489 0.222023 0.335281 0.304257 0.20545 0.184157 0.161903 0.325412 0.300123 0.241809 0.216315 0.329079 0.192159 0.301038 0.290296 0.350127 0.30385 0.33284 0.312083 0.297801 0.237978 0.196205 0.349137 0.232364 0.285143 0.272756 0.302022 0.351807 0.174839 0.302175 0.25245 0.339274 0.222927 0.254692 0.310214 0.282853 0.24827 0.199845 0.297798 0.296549 0.303029 0.328557 0.241543 0.314661 0.385725 -0.027581 0.242124 0.352953 0.238223 0.254982 0.354004 0.338205 0.284555 0.315434 0.30708 0.266865 0.298851 0.336056 0.363471 0.3701 0.16741 0.321012 0.295037 0.417204 0.289055 0.279701 0.267584 0.288166 0.248898 0.2147 0.315381 0.181575 0.293599 0.217744 0.24959 0.276628 0.252833 0.286231 0.259517 0.328284 0.273193 0.272537 0.417374 0.329056 0.233752 0.222025 0.314894 0.294599 0.211834 0.172645 0.240714 0.309123 0.254375 0.276848 0.214229 0.301743 0.341069 0.10271 0.310901 0.212512 0.331638 0.337991 0.338753 0.321264 0.21838 0.340656 0.240702 0.34122 0.322155 0.300499 0.323974 0.217493 0.201413 0.303022 0.196658 0.342765 0.287702 0.242546 0.231014 0.338139 0.327241 0.253901 0.223134 0.322467 0.25106 0.224806 0.224854 0.253981 0.273026 0.198669 0.172635 0.255335 0.288318 0.139313 0.35891 0.27842 0.327671 0.292307 0.244433 0.417109 0.309434 0.337312 0.188259 0.280636 0.311436 0.33408 0.332695 0.272943 0.284644 0.080197 0.260831 0.356373 0.268585 0.364231 0.29018 0.276383 0.258135 0.244105 0.32191 0.320266 0.271522 0.191606 0.279432 0.264892 0.336175 0.327734 0.353917 0.248879 0.301487 0.258516 0.311539 0.327175 0.12896 0.335536 0.350269 0.22583 0.301247 0.190738 0.312494 0.25269 0.275386 0.124118 0.243579 0.330642 0.289719 0.30771 0.325055 0.286336 0.249761 0.299138 0.270883 0.31448 0.243856 0.278602 0.281554 0.329755 0.158548 0.253877 0.314589 0.360739 0.342997 0.27977 0.308291 0.3385 0.284367 0.162514 0.25951 0.29802 0.314572 0.248838 0.283612 0.220921 0.397343 0.313057 0.336191 0.240759 0.238697 0.280036 0.142953 0.186468 0.248215 0.224469 0.310069 0.018946 0.335569 0.313916 0.335718 0.229935 0.269265 0.314046 0.30891 0.333867 0.260068 0.248206 0.313596 0.278258 0.247222 0.358385 0.277176 0.323761 0.067112 0.293972 0.231845 0.259753 0.077582 0.321516 0.270557 0.368308 0.413097 0.252558 0.2456 0.304149 0.318248 0.329794 0.295364 0.346048 0.294515 0.157961 0.348906 0.299827 0.363952 0.100829 0.3135 0.306317 0.071129 0.298442 0.267444 0.324127 0.155831 0.328426 0.339455 0.293682 0.283726 0.261722 0.184087 0.305708 0.256068 0.298811 0.217438 0.265878 0.273605 0.289796 0.320257 0.255085 0.359878 0.302489 0.367869 0.344635 0.324537 0.278985 0.296853 0.349676 0.256949 0.326906 0.368348 0.295257 0.323666 0.230779 0.26031 0.275609 0.264154 0.360468 0.26943 0.352591 0.259197 0.225126 0.30298 0.336715 0.355072 0.368912 0.376919 0.30613 0.26874 0.259427 0.272939 0.072347 0.276096 0.210593 0.30163 0.13782 0.366517 0.293904 0.297652 0.137315 0.264489 0.313308 0.333675 0.351504 0.244889 0.276921 0.273098 0.24445 0.164854 0.363758 0.441923 0.29962 0.302749 0.351003 0.198867 0.317923 0.158193 0.364001 0.187142 0.274109 0.310452 0.218276 0.325068 0.248713 0.333482 0.318803 0.285041 0.290872 0.257508 0.320395 0.320602 0.355112 0.31205 0.31447 0.27943 0.384988 0.302917 0.246191 0.334758 0.389045 0.243338 0.250406 0.275212 0.325911 0.377185 0.308795 0.344669 0.245135 0.29822 0.28973 0.292916 0.280968 0.366705 0.359255 0.316449 0.160962 0.312734 0.205122 0.224652 0.192745 0.223493 0.202811 0.366239 0.263117 0.265179 0.293806 0.287503 0.229977 0.380224 0.356274 0.298063 0.269431 0.117414 0.254 0.245866 0.293257 0.210318 0.13567 0.287545 0.256713 0.287793 0.260519 0.21282 0.195781 0.290501 0.276669 0.320164 0.29927 0.134366 0.30412 0.307737 0.294977 0.340559 0.337168 0.339983 0.310113 0.314679 0.180252 0.27027 0.259754 0.306188 0.093246 0.217766 0.090401 0.385317 0.305922 0.198529 0.337237 0.313863 0.307724 0.219288 0.347587 0.306595 0.282949 0.237665 0.320018 0.315058 0.281436 0.297038 0.354141 0.301376 0.24054 0.279958 0.296369 0.255904 0.295913 0.213225 0.258431 0.284314 0.284916 0.298719 0.232188 0.414336 0.243901 0.157497 0.261809 0.278432 0.330239 0.252029 0.328263 0.313904 0.131408 0.318379 0.299329 0.305537 0.371261 0.321613 0.353756 0.267044 0.270541 0.320642 0.20209 0.257374 0.272811 0.265076 0.166267 0.260762 0.357876 0.291924 0.382356 0.263109 0.13254 0.317241 0.346875 0.217297 0.333787 0.282001 0.268591 0.233263 0.299395 0.237972 0.307068 0.254529 0.243857 0.283326 0.297569 0.216912 0.323068 0.258304 0.24912 0.300801 0.258908 0.260568 0.191644 0.177938 0.323509 0.290291 0.28603 0.207445 0.307908 0.285213 0.26704 0.308182 0.37584 0.235433 0.245345 0.286085 0.300689 0.331115 0.194999 0.296289 0.279823 0.256122 0.132741 0.323538 0.349163 0.197667 0.25934 0.258099 0.190288 0.239842 0.276197 0.289927 0.304352 0.133975 0.317654 0.299708 0.327373 0.280555 0.298098 0.300987 0.378651 0.291815 0.266044 0.250968 0.28556 0.270614 0.280224 0.295499 0.24475 0.300748 0.333824 0.343674 0.235521 0.215244 0.337488 0.263382 0.284794 0.307613 0.320532 0.351376 0.148766 0.30899 0.347251 0.355637 0.301769 0.3331 0.209975 0.253005 0.317373 0.248341 0.265792 0.381929 0.162268 0.310776 0.261983 0.389234 0.201463 0.183376 0.260401 0.30095 0.339398 0.235675 0.339817 0.244362 0.252842 0.263933 0.344553 0.169725 0.362265 0.23349 0.265135 0.368109 0.298368 0.145929 0.238017 0.227869 0.040691 0.350629 0.314775 0.352513 0.328375 0.347639 0.253205 0.30334 0.254653 0.429899 0.363084 0.313349 0.296335 0.244546 0.326767 0.28341 0.210234 0.308868 0.289447 0.300961 0.235914 0.197078 0.380664 0.084 0.25809 0.280802 0.243642 0.256801 0.28812 0.353596 0.296861 0.255499 0.396073 0.247249 0.254268 0.27923 0.295724 0.24084 0.326935 0.353371 0.183517 0.376858 0.273875 0.33112 0.172236 0.331651 0.307431 0.309849 0.29933 0.276056 0.259556 0.303302 0.308613 0.259654 0.329147 0.311691 0.325628 0.241061 0.22522 0.202822 0.273108 0.169137 0.277455 0.219759 0.171576 0.273798 0.311146 0.266051 0.222068 0.289493 0.251724 0.379134 0.265671 0.28352 0.303208 0.193816 0.241029 0.245663 0.362374 0.242447 0.321625 0.374984 0.425967 0.28487 0.201023 0.318578 0.203444 0.376165 0.231353 0.233473 0.253197 0.281687 0.295912 0.353715 0.220409 0.337899 0.240258 0.34488 0.311805 0.265161 0.322775 0.364836 0.307722 0.188934 0.307491 0.286911 0.411284 0.299666 0.329762 0.398818 0.308494 0.253411 0.339377 0.22596 0.23303 0.296089 0.274739 0.242242 0.16914 0.051152 0.210887 0.279617 0.363138 0.25449 0.263888 0.28329 0.330492 0.336827 0.301969 0.221378 0.251001 0.169313 0.337767 0.310161 0.273337 0.236566 0.231717 0.306796 0.417142 0.173959 0.208526 0.313346 0.30219 0.381563 0.305031 0.462976 0.244378 0.346039 0.177017 0.424412 0.326528 0.348065 0.258418 0.373973 0.339991 0.243659 0.215356 0.316093 0.336358 0.332827 0.299602 0.341644 0.319267 0.388058 0.319572 0.298974 0.247563 0.350916 0.298237 0.285188 0.232235 0.329304 0.325293 0.265202 0.218783 0.177851 0.313586 0.236037 0.169047 0.245628 0.302653 0.171973 0.3196 0.261473 0.346561 0.163474 0.378914 0.381615 0.302907 0.337641 0.3073 0.297273 0.126845 0.197773 0.170502 0.292023 0.325708 0.387256 0.327942 0.24439 0.149329 0.302903 0.303385 0.285726 0.300597 0.254075 0.347091 0.212074 0.274742 0.257404 0.272004 0.205699 0.344287 0.293008 0.380494 0.304082 0.342628 0.237985 0.268047 0.336685 0.303682 0.271501 0.278257 0.319727 0.254302 0.315839 0.270584 0.315067 0.307011 0.331533 0.331805 0.307958 0.238322 0.257966 0.420173 0.329225 0.276806 0.238259 0.315989 0.327248 0.271496 0.27063 0.357257 0.293778 0.311293 0.261565 0.380896 0.212724 0.360029 0.311031 0.248078 0.291389 0.283747 0.310662 0.331112 0.339269 0.37813 0.336816 0.374036 0.336993 0.247197 0.275852 0.228873 0.284159 0.347544 0.299015 0.276816 0.185353 0.328335 0.33577 0.230935 0.377006 0.224508 0.135116 0.131978 0.333064 0.359207 0.187918 0.222932 0.242417 0.334745 0.34619 0.410445 0.358476 0.291847 0.222563 0.339049 0.179988 0.194619 0.241703 0.275278 0.321237 0.292544 0.329346 0.196827 0.152711 0.318925 0.347097 0.241007 0.306223 0.317651 0.166939 0.262468 0.281762 0.279141 0.318639 0.290002 0.264813 0.255637 0.258308 0.269466 0.314215 0.278868 0.244804 0.266948 0.242589 0.320569 0.184902 0.283942 0.291043 0.246124 0.264341 0.208154 0.252958 0.214461 0.239911 0.245593 0.239008 0.293333 0.320889 0.353539 0.349063 0.301496 0.208861 0.335605 0.333133 0.216648 0.265162 0.301014 0.224585 0.357315 0.26125 0.275649 0.369829 0.328892 0.320943 0.228957 0.228576 0.207178 0.323844 0.276317 0.274144 0.250503 0.23951 0.30685 0.298015 0.305626 0.209363 0.277775 0.226749 0.304413 0.423024 0.27016 0.255984 0.303753 0.252494 0.18431 0.275167 0.289654 0.22072 0.30121 0.22448 0.378035 0.31653 0.274214 0.172282 0.266842 0.231407 0.287779 0.288378 0.322581 0.332355 0.257144 0.300764 0.207014 0.317011 0.240788 0.2852 0.321864 0.267515 0.307696 0.368644 0.33403 0.24677 0.194236 0.325586 0.228668 0.290511 0.383044 0.239413 0.287795 0.373478 0.352012 0.318758 0.317981 0.320725 0.398167 0.303725 0.295597 0.016663 0.3082 0.252224 0.295819 0.340188 0.272155 0.358964 0.328741 0.259486 0.333821 0.182012 0.277333 0.242602 0.263925 0.266739 0.315308 0.281665 0.334865 0.310459 0.299527 0.347705 0.354573 0.372174 0.320802 0.042284 0.290617 0.232461 0.348166 0.243462 0.301131 0.301524 0.299172 0.364706 0.335491 0.245637 0.293984 0.350276 0.278617 0.3525 0.317295 0.314316 0.318176 0.251056 0.242222 0.321842 0.225792 0.299237 0.260125 0.296436 0.262853 0.323236 0.348211 0.297637 0.274792 0.231894 0.200269 0.273767 0.307572 0.326784 0.308362 0.150609 0.182585 0.311407 0.377249 0.212545 0.265577 0.432698 0.353439 0.328878 0.281438 0.322025 0.319914 0.309736 0.31755 0.322021 0.214672 0.308473 0.317034 0.368319 0.37311 0.285773 0.367952 0.317311 0.277866 0.316163 0.345144 0.225987 0.284056 0.278016 0.340355 0.375798 0.175605 0.16714 0.248146 0.373991 0.258395 0.323975 0.324125 0.227502 0.358681 0.279599 0.370069 0.277815 0.366468 0.191433 0.151653 0.25535 0.202082 0.261977 0.31794 0.313826 0.35267 0.275291 0.257429 0.191414 0.352342 0.309033 0.265969 0.247728 0.304355 0.188222 0.14653 0.415481 0.294915 0.296719 0.256108 0.293024 0.31086 0.248329 0.23727 0.374521 0.294437 0.19607 0.419832 0.280074 0.319973 0.258414 0.3295 0.203595 0.330723 0.352957 0.348269 0.329013 0.242002 0.280438 0.221966 0.339691 0.269209 0.31011 0.285892 0.310783 0.28747 0.202353 0.316519 0.332392 0.280079 0.241468 0.281619 0.2961 0.31016 0.238461 0.244324 0.253134 0.278781 0.29899 0.348357 0.285732 0.294586 0.354097 0.216427 0.340051 0.235705 0.289141 0.093816 0.329966 0.169551 0.287665 0.354 0.294313 0.169419 0.289966 0.254941 0.25414 0.333396 0.263133 0.292036 0.230615 0.365468 0.385395 0.265084 0.251571 0.218594 0.224301 0.379128 0.216551 0.407423 0.326044 0.30923 0.272191 0.190624 0.270358 0.29281 0.140303 0.308359 0.197656 0.172689 0.264773 0.210412 0.054395 0.301511 0.381787 0.307202 0.232968 0.393652 0.28002 0.210257 0.269938 0.329526 0.299106 0.354451 0.279762 0.269252 0.248431 0.271556 0.338602 0.2867 0.318488 0.187306 0.310862 0.266176 0.298368 0.334017 0.298508 0.207763 0.237521 0.358789 0.286361 0.3454 0.266236 0.295085 0.245986 0.286335 0.328451 0.287704 0.298664 0.306651 0.279663 0.301874 0.284896 0.344935 0.349118 0.283614 0.265138 0.150306 0.285047 0.369882 0.277698 0.286658 0.352519 0.238075 0.287266 0.274853 0.149199 0.335081 0.224038 0.283062 0.295301 0.258863 0.291087 0.31188 0.176337 0.268981 0.306054 0.276107 0.19639 0.346958 0.276833 0.352966 0.337843 0.335345 0.252248 0.372415 0.249814 0.262276 0.376847 0.277674 0.310892 0.311001 0.122043 0.310615 0.217692 0.287088 0.328682 0.218028 0.193207 0.263235 0.289787 0.301293 0.27539 0.176984 0.351809 0.261245 0.266308 0.310551 0.314364 0.351356 0.31977 0.291047 0.243183 0.30577 0.262283 0.228374 0.307512 0.245109 0.208966 0.308902 0.246302 0.085709 0.366522 0.251821 0.285031 0.327697 0.389408 0.273301 0.157219 0.25938 0.30315 0.202039 0.345797 0.312365 0.274371 0.334301 0.347494 0.213673 0.331814 0.270049 0.28992 0.300729 0.324254 0.290315 0.293897 0.315887 0.413932 0.164933 0.333539 0.274294 0.375559 0.218775 0.240969 0.234846 0.345211 0.191579 0.300295 0.273552 0.197529 0.284438 0.256902 0.327745 0.237391 0.286582 0.188047 0.183333 0.298344 0.343958 0.311168 0.278453 0.064299 0.297343 0.251006 0.305083 0.295939 0.180681 0.27519 0.322215 0.223208 0.27913 0.368439 0.185491 0.340768 0.233591 0.306209 0.281357 0.273818 0.371939 0.30056 0.210663 0.335025 0.292739 0.209374 0.140516 0.325313 0.273262 0.396008 0.227918 0.335841 0.23007 0.272149 0.326864 0.281893 0.271796 0.308734 0.343696 0.33988 0.201991 0.32092 0.276606 0.286667 0.320079 0.25158 0.305025 0.289363 0.266834 0.306152 0.267289 0.25805 0.297614 0.340092 0.272187 0.393732 0.292525 0.289037 0.332047 0.219501 0.232956 0.324246 0.276142 0.354401 0.377863 0.353466 0.399487 0.142949 0.16202 0.230108 0.385134 0.28833 0.284697 0.20053 0.18276 0.244487 0.342043 0.142439 0.267242 0.258628 0.061097 0.391374 0.288457 0.361091 0.409928 0.310675 0.295404 0.193796 0.324526 0.295892 0.317491 0.267753 0.137185 0.289676 0.346445 0.287917 0.227232 0.276049 0.311959 0.164824 0.276782 0.243274 0.208932 0.325096 0.359694 0.258932 0.345872 0.267419 0.294106 0.265551 0.283505 0.321878 0.278597 0.277737 0.250989 0.296418 0.260572 0.292122 0.296721 0.276513 0.271778 0.260606 0.358328 0.249477 0.20359 0.263529 0.210108 0.278635 0.246271 0.261836 0.333228 0.180907 0.296859 0.301962 0.17483 0.27891 0.145101 0.289476 0.359431 0.33135 0.305624 0.353597 0.105161 0.327878 0.155886 0.302805 0.357971 0.183768 0.267101 0.286475 0.261638 0.150678 0.160231 0.341274 0.313416 0.11573 0.185373 0.339597 0.284648 0.323198 0.3319 0.214937 0.223981 0.310372 0.317931 0.355718 0.224933 0.35927 0.331852 0.269139 0.279592 0.338739 0.341238 0.341765 0.291575 0.288987 0.256552 0.274548 0.24295 0.287558 0.223294 0.317058 0.368656 0.271857 0.211007 0.265195 0.24965 0.249967 0.238827 0.209148 0.322111 0.267338 0.372302 0.263771 0.161937 0.232052 0.281354 0.403266 0.297435 0.121358 0.187151 0.330471 0.287572 0.179709 0.216632 0.309756 0.277536 0.2237 0.228454 0.191964 0.269035 0.358954 0.292891 0.284444 0.311148 0.319556 0.299106 0.249746 0.269821 0.32963 0.240157 0.378265 0.30704 0.342775 0.151293 0.275717 0.290567 0.227038 0.365877 0.333961 0.276769 0.337041 0.395267 0.325123 0.208007 0.300194 0.319773 0.309863 0.265784 0.258628 0.20191 0.139022 0.268043 0.31568 0.203185 0.315158 0.346817 0.211017 0.304134 0.113446 0.269804 0.26048 0.31393 0.285656 0.269039 0.197965 0.306045 0.346572 0.280313 0.336102 0.202394 0.290503 0.199719 0.29992 0.331596 0.171181 0.319364 0.314239 0.294785 0.342237 0.30895 0.262242 0.321522 0.311722 0.397746 0.32576 0.319008 0.325617 0.267249 0.299574 0.295706 0.174904 0.163125 0.136295 0.294018 0.352364 0.224941 0.229878 0.315693 0.244746 0.201405 0.339572 0.322319 0.301522 0.248495 0.234305 0.379356 0.2599 0.309366 0.138336 0.217198 0.216322 0.306934 0.251066 0.300171 0.289127 0.305697 0.260782 0.221968 0.205098 0.279179 0.39221 0.347438 0.248883 0.122581 0.266427 0.326346 0.155456 0.258428 0.31939 0.204307 0.326518 0.381621 0.250018 0.189551 0.31932 0.350327 0.316617 0.267578 0.344201 0.275629 0.262239 0.300355 0.256949 0.327508 0.319124 0.284633 0.227291 0.306744 0.231581 0.479004 0.225334 0.293949 0.229175 0.321628 0.284721 0.325135 0.248598 0.285491 0.21687 0.262479 0.285838 0.249268 0.353277 0.34064 0.323863 0.22527 0.265954 0.271282 0.356396 0.257863 0.380584 0.287385 0.233036 0.207952 0.287131 0.357151 0.310568 0.277627 0.334784 0.140234 0.319592 0.245267 0.130971 0.231217 0.272138 0.199429 0.240957 0.362268 0.25241 0.254355 0.280761 0.274116 0.235249 0.200305 0.330725 0.189741 0.387007 0.19343 0.252122 0.366403 0.281413 0.280628 0.327483 0.170318 0.327983 0.303599 0.299808 0.280995 0.34053 0.285312 0.259248 0.215491 0.263782 0.292672 0.248216 0.22249 0.264344 0.151535 0.345298 0.294617 0.323974 0.249915 0.291923 0.26132 0.340145 0.156106 0.316455 0.224818 0.318829 0.330572 0.237218 0.275613 0.264615 0.316323 0.282946 0.263649 0.291578 0.210894 0.274263 0.325221 0.252856 0.255134 0.297845 0.238956 0.341969 0.208966 0.286386 0.310598 0.357019 0.259844 0.253034 0.190334 0.184634 0.254194 0.272003 0.329586 0.289793 0.220541 0.309195 0.340271 0.280205 0.200144 0.201839 0.37314 0.290148 0.315126 0.363831 0.239742 0.285868 0.221932 0.177001 0.343621 0.304421 0.235243 0.375249 0.322256 0.279728 0.283326 0.321387 0.306968 0.331412 0.227597 0.374951 0.327291 0.331528 0.28547 0.265995 0.309083 0.290093 0.297622 0.314989 0.343249 0.327576 0.260444 0.249814 0.312663 0.377113 0.307264 0.229414 0.318848 0.277259 0.31789 0.322813 0.363772 0.299392 0.307298 0.319606 0.236833 0.344557 0.071317 0.23473 0.298199 0.196593 0.275641 0.267768 0.283771 0.238425 0.341821 0.251917 0.323777 0.308906 0.261449 0.264241 0.287329 0.307109 0.379684 0.203935 0.357644 0.274652 0.294922 0.310168 0.244806 0.284772 0.246636 0.355057 0.296417 0.237876 0.304142 0.295944 0.346964 0.309519 0.322843 0.269524 0.319778 0.310266 0.248717 0.270815 0.239658 0.277056 0.314774 0.313108 0.286895 0.298219 0.274563 0.229262 0.324873 0.304559 0.238834 0.15788 0.311606 0.136881 0.36402 0.265867 0.255552 0.267808 0.381064 0.254365 0.279545 0.337648 0.256503 0.25335 0.28767 0.242781 0.228449 0.364669 0.402384 0.29111 0.275911 0.386377 0.306527 0.232978 0.29013 0.294552 0.179748 0.246614 0.344225 0.092403 0.392433 0.317003 0.329143 0.141273 0.349156 0.272885 0.251513 0.300049 0.358175 0.301881 0.326637 0.394309 0.220702 0.25024 0.347021 0.400447 0.280106 0.199254 0.299357 0.29808 0.311875 0.232943 0.403662 0.227287 0.234406 0.338429 0.276033 0.28996 0.294313 0.25021 0.248514 0.25131 0.297311 0.346233 0.266824 0.304218 0.40583 0.303457 0.250592 0.293654 0.300624 0.285043 0.421678 0.206921 0.277495 0.338462 0.352866 0.04786 0.240137 0.376985 0.230767 0.264265 0.501285 0.298917 0.30444 0.292983 0.277387 0.303061 0.231764 0.371141 0.285824 0.318133 0.282675 0.270856 0.324402 0.261163 0.390534 0.19934 0.234427 0.288679 0.351825 0.292555 0.280275 0.366754 0.242378 0.281961 0.293931 0.205903 0.295189 0.213902 0.281183 0.291573 0.236562 0.22093 0.329887 0.192958 0.300672 0.131114 0.217739 0.28414 0.327594 0.194091 0.293265 0.317357 0.350055 0.364437 0.309686 0.269504 0.272122 0.317219 0.328001 0.159644 0.331514 0.252321 0.308245 0.271266 0.270144 0.311314 0.224716 0.369109 0.385464 0.107169 0.090704 0.348093 0.292029 0.323992 0.28038 0.307389 0.358111 0.283041 0.353466 0.192209 0.284526 0.322656 0.308839 0.376929 0.402423 0.341127 0.188334 0.2813 0.208547 0.231263 0.215065 0.347869 0.259046 0.306971 0.321308 0.264266 0.222634 0.215925 0.240314 0.196566 0.311771 0.284952 0.383778 0.309069 0.28483 0.358813 0.289131 0.310464 0.357121 0.242647 0.368257 0.340026 0.226065 0.3465 0.150752 0.271348 0.31391 0.211077 0.202254 0.300874 0.238063 0.325978 0.308963 0.262052 0.216742 0.391455 0.344009 0.44179 0.072184 0.293973 0.414652 0.296657 0.266627 0.360569 0.44052 0.208657 0.319061 0.262123 0.347528 0.292161 0.288317 0.259294 0.327319 0.292185 0.267816 0.283843 0.369966 0.28754 0.306078 0.335084 0.335419 0.266584 0.304182 0.228143 0.344272 0.338399 0.337052 0.277192 0.282388 0.341719 0.324216 0.340253 0.340478 0.17741 0.280516 0.259917 0.228813 0.26806 0.320781 0.241275 0.297246 0.379107 0.347458 0.273146 0.38265 0.346407 0.309039 0.284317 0.196433 0.30596 0.283678 0.290629 0.297022 0.176376 0.291129 0.307044 0.198183 0.343655 0.345401 0.291216 0.243802 0.312108 0.307884 0.28611 0.288601 0.21245 0.217744 0.324415 0.270204 0.227413 0.305325 0.272748 0.330979 0.315865 0.317568 0.240293 0.301411 0.286851 0.241438 0.226587 0.310494 0.324071 0.189669 0.300596 0.329341 0.233103 0.313184 0.314475 0.284453 0.107931 0.301711 0.263338 0.329102 0.321859 0.242878 0.223662 0.274164 0.313924 0.274456 0.297351 0.290954 0.224721 0.278797 0.331562 0.326167 0.30322 0.279267 0.207409 0.356632 0.310873 0.238774 0.31086 0.249085 0.276965 0.323949 0.326713 0.372422 0.292265 0.383475 0.255336 0.282281 0.280039 0.201469 0.35382 0.280808 0.330027 0.342164 0.330232 0.060757 0.307661 0.280219 0.2921 0.309622 0.215592 0.356588 0.333358 0.319597 0.266863 0.257576 0.286534 0.382872 0.226559 0.267322 0.15878 0.254204 0.307498 0.27457 0.338207 0.234774 0.267049 0.225753 0.332602 0.106009 0.153357 0.143318 0.409229 0.236156 0.306514 0.322831 0.303948 0.109855 0.337152 0.108805 0.265098 0.258038 0.322868 0.334828 0.262824 0.372448 0.04476 0.287751 0.373474 0.201265 0.322811 0.320048 0.319921 0.313438 0.310287 0.304758 0.212291 0.30892 0.241444 0.358088 0.140545 0.327986 0.261798 0.24155 0.333741 0.261534 0.291478 0.299857 0.30472 0.295565 0.220304 0.341478 0.24337 0.156608 0.225785 0.391344 0.281072 0.262269 0.375625 0.318758 0.291957 0.297737 0.257106 0.403935 0.308783 0.2156 0.264317 0.316905 0.343725 0.241726 0.328938 0.275896 0.291014 0.315422 0.34641 0.281359 0.399011 0.353153 0.270007 0.264306 0.20069 0.286786 0.301307 0.240933 0.361882 0.285435 0.229057 0.196216 0.325429 0.209848 0.311303 0.279401 0.308471 0.355384 0.236462 0.242149 0.315694 0.122678 0.260376 0.280328 0.237028 0.35067 0.213556 0.321279 0.226462 0.169583 0.459772 0.419341 0.303486 0.312984 0.33364 0.264667 0.322089 0.325275 0.246773 0.332721 0.30809 0.265695 0.299838 0.350284 0.255734 0.291757 0.278032 0.291684 0.253957 0.278756 0.307738 0.299556 0.15845 0.31367 0.357487 0.303373 0.23653 0.251336 0.222971 0.298175 0.390724 0.29145 0.30051 0.313651 0.308814 0.27339 0.329489 0.265822 0.220991 0.248721 0.297976 0.265239 0.238248 0.29285 0.179766 0.225736 0.197523 0.289145 0.189753 0.307621 0.215435 0.297186 0.330276 0.297456 0.292662 0.281948 0.369841 0.233268 0.327089 0.386139 0.28111 0.232728 0.369164 0.316161 0.341687 0.215098 0.280949 0.227552 0.249014 0.2267 0.196511 0.259528 0.303434 0.1715 0.261038 0.26672 0.35946 0.28058 0.336688 0.308455 0.20564 0.275655 0.256154 0.25273 0.066242 0.298624 0.352398 0.193957 0.317189 0.24647 0.248123 0.260735 0.273432 0.27345 0.35688 0.326615 0.322328 0.290438 0.305846 0.283618 0.308833 0.304975 0.303907 0.295469 0.307011 0.335936 0.269951 0.146432 0.288413 0.306998 0.274497 0.294069 0.214961 0.285695 0.192275 0.285914 0.313208 0.365728 0.328627 0.259747 0.36127 0.325405 0.271002 0.26347 0.324391 0.278602 0.203722 0.291883 0.273395 0.296723 0.218365 0.268699 0.299881 0.348146 0.256353 0.206754 0.251521 0.238575 0.371188 0.319432 0.432666 0.319114 0.321151 0.198967 0.330092 0.127683 0.316847 0.387594 0.190644 0.339556 0.360585 0.35322 0.336269 0.287268 0.269698 0.242128 0.287257 0.291707 0.300514 0.376466 0.321568 0.278868 0.32055 0.266164 0.307728 0.288388 0.243098 0.306786 0.291566 0.300896 0.145172 0.158347 0.178359 0.307938 0.377615 0.190986 0.268903 0.237211 0.289012 0.250773 0.289669 0.227459 0.249686 0.230355 0.378271 0.191374 0.274981 0.349219 0.271935 0.314425 0.256256 0.259038 0.306934 0.360933 0.201471 0.259538 0.248973 0.244175 0.326141 0.317375 0.240854 0.314648 0.175155 0.207365 0.392535 0.251986 0.336564 0.344139 0.316717 0.323405 0.335135 0.317244 0.375845 0.183107 0.333785 0.297166 0.23972 0.199607 0.348365 0.326404 0.32386 0.343013 0.304184 0.276436 0.265833 0.268272 0.318495 0.159643 0.329215 0.300997 0.271865 0.307969 0.278302 0.114202 0.252887 0.132272 0.259024 0.355017 0.350561 0.244479 0.32187 0.371124 0.290157 0.332662 0.314505 0.296123 0.314522 0.262896 0.305844 0.265819 0.179463 0.382181 0.302438 0.404882 0.337264 0.332987 0.328044 0.324035 0.278702 0.252959 0.307192 0.197456 0.252686 0.316249 0.298322 0.340144 0.295945 0.2465 0.267066 0.255819 0.262511 0.110337 0.209949 0.26438 0.004085 0.336601 0.323775 0.262516 0.307304 0.317924 0.261488 0.388976 0.201404 0.316762 0.261911 0.304985 0.313279 0.290405 0.321435 0.335064 0.317263 0.277485 0.304745 0.290578 0.32692 0.301708 0.172651 0.306671 0.329637 0.29542 0.256401 0.196146 0.167483 0.312257 0.303516 0.243545 0.242082 0.263977 0.280811 0.282345 0.175505 0.436189 0.347098 0.24307 0.334071 0.298558 0.197592 0.262447 0.29961 0.352765 0.325999 0.301061 0.328134 0.327457 0.270992 0.302153 0.323795 0.259927 0.216059 0.268682 0.286841 0.255806 0.288718 0.318353 0.316098 0.308402 0.263933 0.22092 0.244433 0.314913 0.347747 0.293365 0.291326 0.218348 0.317699 0.262564 0.325364 0.23854 0.27305 0.29451 0.347924 0.2721 0.137059 0.294671 0.276574 0.251143 0.284496 0.288089 0.237776 0.118743 0.288684 0.330655 0.114273 0.275419 0.344748 0.264979 0.350199 0.310142 0.308678 0.436808 0.151382 0.273288 0.207663 0.32603 0.377625 0.250304 0.333598 0.234195 0.279205 0.319255 0.26992 0.300773 0.310412 0.266969 0.273364 0.351307 0.176247 0.303536 0.173614 0.173171 0.18521 0.298921 0.249669 0.249471 0.35154 0.302889 0.268912 0.195098 0.177973 0.276469 0.238699 0.265658 0.315543 0.302223 0.306542 0.321831 0.270164 0.278948 0.362502 0.316933 0.301594 0.199717 0.191793 0.273436 0.343274 0.262114 0.2893 0.333131 0.262736 0.248746 0.296928 0.270429 0.222476 0.280129 0.240655 0.365455 0.29671 0.295542 0.145308 0.243464 0.211316 0.265139 0.118337 0.300876 0.337617 0.309043 0.182801 0.330382 0.317625 0.358521 0.311522 0.305592 0.317343 0.329078 0.310719 0.278398 0.169844 0.236269 0.281842 0.320388 0.31125 0.253535 0.292786 0.249908 0.332443 0.313166 0.35156 0.291017 0.332265 0.170916 0.294552 0.445771 0.26628 0.230423 0.321001 0.304009 0.258403 0.232544 0.098292 0.300187 0.186823 0.201544 0.252675 0.158077 0.281384 0.237247 0.229955 0.247138 0.26668 0.35525 0.281182 0.316746 0.354398 0.274291 0.36398 0.344389 0.323588 0.260257 0.303037 0.289524 0.289399 0.320829 0.312848 0.302171 0.31513 0.231961 0.261017 0.34294 0.354426 0.310179 0.208415 0.295005 0.270816 0.291333 0.340533 0.300698 0.262333 0.319575 0.31517 0.356613 0.382027 0.317796 0.104486 0.377177 0.201371 0.261897 0.249919 0.34795 0.308616 0.305812 0.27908 0.305461 0.29422 0.234275 0.298442 0.326842 0.372355 0.348152 0.342574 0.188317 0.351165 0.341144 0.223416 0.337676 0.364538 0.279231 0.248763 0.372755 0.307024 0.295217 0.270917 0.276185 0.332734 0.274851 0.275582 0.289849 0.362521 0.371327 0.260808 0.288656 0.312287 0.347178 0.247143 0.211566 0.25578 0.204558 0.294091 0.34139 0.347572 0.300973 0.308138 0.293621 0.331193 0.268398 0.31027 0.18297 0.29134 0.317305 0.333631 0.437463 0.337559 0.291544 0.240445 0.278646 0.281696 0.3843 0.263158 0.263514 0.330396 0.333842 0.243128 0.245173 0.273671 0.295809 0.33704 0.37001 0.25423 0.303685 0.295968 0.135822 0.255978 0.313014 0.275513 0.264811 0.296781 0.267196 0.304876 0.340475 0.03944 0.228837 0.309268 0.407666 0.214475 0.18106 0.334313 0.284532 0.317736 0.333181 0.315256 0.198853 0.307437 0.227099 0.410776 0.172915 0.233889 0.357879 0.260013 0.257053 0.321865 0.329839 0.357479 0.318807 0.331464 0.346664 0.319539 0.304802 0.308472 0.28987 0.312182 0.178162 0.246627 0.271321 0.270036 0.199017 0.346165 0.230172 0.267598 0.331579 0.272838 0.241235 0.23379 0.270897 0.375006 0.432858 0.136541 0.275376 0.353241 0.261682 0.352032 0.289505 0.126076 0.207656 0.221014 0.250526 0.33676 0.209899 0.382431 0.312297 0.335715 0.299949 0.31971 0.274785 0.268771 0.247621 0.359022 0.379787 0.227228 0.281447 0.333402 0.294184 0.237059 0.289866 0.305684 0.24591 0.322087 0.207047 0.283689 0.330653 0.325744 0.302931 0.323907 0.260651 0.321056 0.241278 0.203465 0.186079 0.298353 0.359445 0.283562 0.222825 0.178694 0.255313 0.253824 0.292064 0.226655 0.284805 0.324546 0.247802 0.214082 0.339125 0.289609 0.178726 0.264192 0.348285 0.358023 0.327311 0.276863 0.297321 0.262935 0.234393 0.250604 0.379041 0.31503 0.270191 0.250239 0.336823 0.357802 0.320826 0.257167 0.293941 0.285097 0.319571 0.251747 0.292279 0.260223 0.255438 0.24511 0.420639 0.193628 0.300649 0.370828 0.370588 0.317077 0.23335 0.373003 0.355487 0.309837 0.302384 0.309048 0.147564 0.254548 0.323543 0.275133 0.366661 0.274253 0.295089 0.244021 0.270948 0.322081 0.286558 0.323679 0.290078 0.255718 0.301401 0.206314 0.286222 0.228349 0.275575 0.334492 0.348085 0.325095 0.327444 0.345328 0.231616 0.313377 0.341242 0.360075 0.263562 0.304289 0.325929 0.250524 0.366621 0.299399 0.2476 0.216695 0.140629 0.279855 0.36649 0.251264 0.369555 0.297244 0.307397 0.185096 0.302588 0.514318 0.308186 0.329102 0.276194 0.369627 0.269399 0.325977 0.286505 0.311927 0.309696 0.304293 0.12852 0.256238 0.157847 0.26162 0.371005 0.207223 0.309875 0.295737 0.318914 0.193193 0.265129 0.311601 0.279258 0.280042 0.27628 0.247929 0.279602 0.458385 0.286166 0.316757 0.22042 0.280271 0.231771 0.31924 0.230827 0.38918 0.296705 0.268717 0.271427 0.302931 0.285126 0.269318 0.258631 0.421325 0.166285 0.26583 0.333929 0.169792 0.336965 0.242666 0.318607 0.293644 0.348335 0.193274 0.362211 0.214403 0.213859 0.376433 0.256364 0.251796 0.342134 0.234704 0.264039 0.302898 0.243904 0.172602 0.313053 0.260072 0.267733 0.262618 0.316857 0.351104 0.243497 0.26109 0.353674 0.305833 0.178568 0.355765 0.349552 0.211998 0.363787 0.250214 0.363951 0.108956 0.282859 0.304647 0.323803 0.329161 0.324102 0.303993 0.264347 0.289445 0.306092 0.280992 0.253339 0.346212 0.255124 0.323278 0.253034 0.330848 0.358954 0.280879 0.203112 0.377374 0.39686 0.302736 0.202919 0.191627 0.217414 0.27396 0.201306 0.29169 0.310216 0.321083 0.318836 0.274493 0.367802 0.305242 0.255128 0.303478 0.299614 0.26704 0.281336 0.309013 0.361214 0.145789 0.244217 0.331762 0.261825 0.259708 0.299076 0.244813 0.390263 0.203397 0.316936 0.337721 0.25839 0.241727 0.332741 0.343344 0.258038 0.275494 0.240183 0.226342 0.279043 0.362239 0.287464 0.228952 0.283626 0.345165 0.309177 0.197721 0.307081 0.466833 0.267658 0.160804 0.303104 0.313918 0.340201 0.299074 0.266315 0.111879 0.235125 0.320339 0.219366 0.229394 0.305551 0.262847 0.249254 0.32649 0.325428 0.199113 0.297485 0.292468 0.19818 0.367199 0.260348 0.287571 0.295357 0.359583 0.249624 0.208098 0.278476 0.310546 0.292001 0.263642 0.33936 0.278794 0.415595 0.289236 0.365676 0.323574 0.264283 0.326588 0.21852 0.34104 0.346183 0.186114 0.302218 0.310875 0.281909 0.193327 0.300933 0.243842 0.200934 0.326471 0.277911 0.298449 0.319224 0.234519 0.309795 0.178089 0.275685 0.460743 0.250086 0.267522 0.317116 0.290954 0.327815 0.255741 0.355941 0.204617 0.344872 0.274158 0.220305 0.356414 0.273176 0.257696 0.253041 0.295652 0.31639 0.333696 0.318338 0.108276 0.236946 0.299553 0.2406 0.233022 0.292623 0.276432 0.239282 0.269599 0.25932 0.291687 0.18365 0.294263 0.319824 0.320577 0.282223 0.308074 0.390886 0.234187 0.331889 0.307683 0.298685 0.306584 0.295111 0.321404 0.296302 0.178525 0.261195 0.279449 0.210508 0.205256 0.357588 0.276455 0.266135 0.338832 0.28622 0.319496 0.293 0.24874 0.322482 0.248705 0.270641 0.252702 0.293911 0.322489 0.307606 0.324799 0.305149 0.168276 0.244732 0.325842 0.275065 0.336874 0.153845 0.304128 0.285117 0.347244 0.318191 0.286675 0.32551 0.303958 0.324408 0.257275 0.344957 0.244811 0.337817 0.359852 0.223183 0.329201 0.339276 0.315043 0.294068 0.326275 0.308283 0.046724 0.245097 0.212912 0.277362 0.3062 0.263356 0.322651 0.278951 0.285691 0.262122 0.380705 0.238704 0.148139 0.232975 0.330933 0.31008 0.254461 0.289023 0.319965 0.295235 0.32187 0.25501 0.232813 0.318092 0.36337 0.212349 0.291157 0.286172 0.26041 0.214676 0.21079 0.214521 0.322661 0.278973 0.275948 0.260241 0.22861 0.31872 0.326051 0.293905 0.303071 0.272877 0.266985 0.266483 0.35605 0.277967 0.29651 0.225877 0.333754 0.218833 0.374282 0.422055 0.325652 0.257959 0.30299 0.294846 0.289848 0.191263 0.308747 0.244502 0.260309 0.259231 0.273232 0.255912 0.329055 0.192158 0.360042 0.411634 0.307363 0.273239 0.248824 0.38834 0.342532 0.326824 0.210033 0.244122 0.311791 0.31209 0.340335 0.20784 0.257317 0.286422 0.323232 0.21978 0.363096 0.344936 0.339388 0.305167 0.257972 0.274108 0.345955 0.205372 0.313071 0.278043 0.255739 0.232371 0.294871 0.330213 0.340563 0.33594 0.263854 0.302074 0.281928 0.338438 0.274962 0.198382 0.293794 0.241043 0.185412 0.374377 0.30152 0.283213 0.382369 0.360264 0.267707 0.281641 0.287121 0.369722 0.183593 0.349179 0.261284 0.29495 0.287234 0.30317 0.283463 0.296776 0.308806 0.356995 0.353893 0.391433 0.232161 0.320365 0.312447 0.248964 0.215327 0.236614 0.320837 0.290108 0.267507 0.267433 0.195861 0.269413 0.303307 0.235064 0.264983 0.31614 0.389158 0.309743 0.309015 0.241687 0.340032 0.285209 0.291783 0.317877 0.183916 0.275687 0.174929 0.330756 0.267569 0.352979 0.274479 0.312161 0.335169 0.236761 0.239091 0.253713 0.199581 0.295106 0.349585 0.182306 0.249718 0.185531 0.254097 0.255047 0.281077 0.244791 0.300346 0.213259 0.29599 0.35065 0.198616 0.235527 0.266577 0.208707 0.219777 0.211373 0.291459 0.229265 0.128305 0.305279 0.303236 0.255261 0.334563 0.187639 0.388761 0.260049 0.248642 0.147891 0.351793 0.31627 0.127405 0.243398 0.326485 0.30641 0.244335 0.297063 0.229778 0.383524 0.157555 0.228899 0.25426 0.270492 0.321598 0.179033 0.269596 0.326354 0.244185 0.286203 0.316605 0.149154 0.244365 0.329595 0.316806 0.271376 0.297357 0.234914 0.257254 0.392584 0.26664 0.305507 0.352549 0.274099 0.262776 0.300699 0.239412 0.329645 0.31223 0.200573 0.410456 0.319863 0.431726 0.322993 0.311836 0.317268 0.266843 0.319497 0.279983 0.339852 0.25292 0.297751 0.303872 0.342657 0.204467 0.287093 0.178947 0.287882 0.271034 0.3103 0.252724 0.259027 0.23788 0.184726 0.253681 0.229441 0.326215 0.108379 0.306614 0.325815 0.248297 0.317307 0.282076 0.27963 0.212118 0.16839 0.334938 0.337073 0.29216 0.365227 0.367607 0.280301 0.285461 0.294217 0.276079 0.271405 0.292741 0.276038 0.287202 0.284514 0.322059 0.124664 0.381309 0.407567 0.214959 0.245611 0.258852 0.364791 0.269954 0.272598 0.378049 0.195745 0.360887 0.333842 0.199027 0.223557 0.354885 0.37946 0.308837 0.284076 0.256316 0.268014 0.341447 0.347469 0.299749 0.353563 0.279528 0.294104 0.278193 0.284556 0.206408 0.331686 0.304857 0.261398 0.307505 0.320454 0.216988 0.175356 0.232662 0.250175 0.287129 0.339503 0.263355 0.242511 0.297256 0.240137 0.307939 0.35652 0.321558 0.26059 0.264953 0.274469 0.336447 0.190428 0.098621 0.117586 0.291589 0.323531 0.247271 0.366234 0.340393 0.315646 0.245284 0.32074 0.229524 0.337076 0.278246 0.318286 0.316502 0.297908 0.281236 0.388561 0.326943 0.375853 0.356333 0.276803 0.278381 0.285408 0.264663 0.231704 0.26906 0.389249 0.23749 0.350777 0.296614 0.235236 0.230251 0.351584 0.31832 0.32 0.206546 0.391205 0.275905 0.206347 0.268817 0.335661 0.411758 0.289997 0.385992 0.28655 0.288797 0.264844 0.275592 0.185556 0.327756 0.311362 0.270494 0.26476 0.229888 0.342163 0.271439 0.360055 0.295298 0.416384 0.302795 0.332689 0.22294 0.283701 0.313485 0.328001 0.344603 0.159635 0.391566 0.174111 0.322861 0.335798 0.127932 0.152692 0.260291 0.363623 0.414775 0.357786 0.365454 0.353686 0.287249 0.309659 0.211231 0.290159 0.324496 0.326767 0.323475 0.151701 0.26985 0.196185 0.300079 0.406816 0.22282 0.21275 0.284776 0.309272 0.233768 0.294751 0.321761 0.273771 0.285634 0.307328 0.271263 0.091253 0.346582 0.28627 0.219133 0.247729 0.309732 0.257804 0.269369 0.31 0.187482 0.275758 0.352856 0.077212 0.253413 0.307436 0.318825 0.291832 0.366051 0.269287 0.189673 0.177774 0.306258 0.283864 0.314201 0.269066 0.360648 0.097903 0.251767 0.195892 0.287898 0.285453 0.310112 0.220226 0.211459 0.462643 0.124351 0.233328 0.287771 0.320536 0.324323 0.314428 0.230697 0.280473 0.334534 0.293425 0.21309 0.31494 0.326012 0.231167 0.264627 0.295657 0.297848 0.235065 0.285636 0.16859 0.250473 0.297734 0.297531 0.291888 0.304314 0.166321 0.353331 0.284035 0.255241 0.225919 0.35564 0.278845 0.386158 0.246728 0.28457 0.332994 0.28535 0.274649 0.215767 0.208483 0.265528 0.426227 0.25644 0.323805 0.33715 0.259781 0.299125 0.210521 0.325453 0.303886 0.297538 0.316625 0.181172 0.263934 0.299717 0.227794 0.314622 0.316337 0.179742 0.292779 0.10441 0.251225 0.340537 0.414159 0.257816 0.17703 0.320178 0.267199 0.358413 0.337139 0.243072 0.36016 0.239925 0.288377 0.30231 0.216115 0.318477 0.298752 0.23918 0.337933 0.299038 0.262061 0.325664 0.30931 0.246046 0.194482 0.29635 0.260313 0.346584 0.326847 0.349109 0.22026 0.247702 0.359738 0.379281 0.158491 0.259928 0.257039 0.188928 0.26432 0.185957 0.369131 0.239359 0.256554 0.200837 0.30856 0.269015 0.288324 0.258981 0.309849 0.285995 0.339687 0.304156 0.224545 0.336465 0.301085 0.343216 0.381031 0.173346 0.286636 0.282502 0.348559 0.253229 0.368262 0.319313 0.303243 0.271418 0.247249 0.219372 0.312718 0.283489 0.458834 0.324173 0.377011 0.277164 0.297236 0.259633 0.274101 0.225546 0.290007 0.397759 0.297884 0.305818 0.367011 0.38035 0.314216 0.312451 0.275703 0.288186 0.290224 0.189042 0.297331 0.32605 0.253292 0.328832 0.252571 0.223307 0.31543 0.324034 0.254789 0.274353 0.338501 0.204934 0.299156 0.221199 0.22664 0.227972 0.317692 0.236827 0.19811 0.348123 0.32168 0.353013 0.265736 0.158018 0.203295 0.298999 0.363487 0.168466 0.330935 0.238531 0.250632 0.379901 0.257645 0.078056 0.229734 0.285881 0.295933 0.299611 0.274643 0.16025 0.294205 0.318152 0.196438 0.358111 0.317383 0.312205 0.310283 0.316587 0.293632 0.3247 0.26317 0.190317 0.298064 0.265825 0.227715 0.383935 0.320937 0.298566 0.258367 0.328425 0.298021 0.360459 0.213939 0.327075 0.251383 0.429529 0.225634 0.34987 0.333503 0.389192 0.310287 0.401909 0.369532 0.249782 0.364966 0.31504 0.330066 0.113603 0.263679 0.237789 0.15802 0.304368 0.281897 0.244684 0.405792 0.326225 0.352454 0.231166 0.262778 0.196416 0.319786 0.205813 0.256158 0.272995 0.089332 0.213073 0.261048 0.288809 0.368149 0.240493 0.303115 0.452507 0.220015 0.152038 0.225888 0.132382 0.387085 0.348146 0.259855 0.037426 0.280547 0.230113 0.333303 0.271902 0.314932 0.242565 0.023289 0.293753 0.255932 0.259895 0.30708 0.169682 0.288707 0.295014 0.244551 0.153584 0.272834 0.386789 0.278085 0.322672 0.225909 0.284423 0.116256 0.27582 0.347989 0.329082 0.242605 0.31064 0.21554 0.250542 0.302892 0.226834 0.256986 0.359439 0.302225 0.244484 0.333199 0.283743 0.221553 0.254229 0.446237 0.343156 0.253966 0.229201 0.256159 0.346256 0.130313 0.31365 0.230328 0.320698 0.227828 0.340118 0.138159 0.260241 0.312986 0.409592 0.223141 0.307888 0.259123 0.289659 0.320392 0.194209 0.345853 0.342504 0.244061 0.303959 0.27484 0.340801 0.174682 0.316877 0.230444 0.01901 0.278873 0.248309 0.295423 0.315268 0.239702 0.277074 0.11937 0.261362 0.235294 0.295371 0.109538 0.314953 0.04151 0.20767 0.216116 0.279474 0.238819 0.273631 0.179176 0.281951 0.222796 0.168573 0.231829 0.28469 0.306366 0.304653 0.285345 0.280256 0.262803 0.165072 0.274613 0.263169 0.296195 0.217863 0.266564 0.309491 0.298165 0.238305 0.278641 0.06627 0.192428 0.285967 0.339963 0.222636 0.27542 0.325787 0.174073 0.369666 0.178826 0.31436 0.346672 0.263579 0.261417 0.281549 0.155119 0.29902 0.334836 0.285301 0.249447 0.283028 0.341098 0.360705 0.286176 0.301369 0.147843 0.273327 0.288021 0.280233 0.32958 0.416324 0.275821 0.354875 0.341093 0.174549 0.304615 0.30722 0.320926 0.283518 0.267472 0.270658 0.326506 0.308049 0.282053 0.321493 0.312177 0.308487 0.209083 0.342879 0.290416 0.316531 0.30387 0.33727 0.360551 0.313283 0.299356 0.299897 0.28648 0.29167 0.399415 0.24929 0.268708 0.222182 0.362423 0.261981 0.270905 0.360679 0.276782 0.27884 0.259276 0.370951 0.315462 0.251164 0.305122 0.360549 0.318958 0.236692 0.282819 0.297489 0.228399 0.361716 0.301516 0.280158 0.331986 0.261541 0.203997 0.288758 0.381446 0.304373 0.374876 0.307179 0.237259 0.212373 0.214501 0.277648 0.226968 0.329221 0.278963 0.321258 0.344235 0.309799 0.230063 0.260993 0.286621 0.218812 0.299319 0.301086 0.242972 0.331074 0.296141 0.266463 0.021364 0.272584 0.386892 0.358969 0.179553 0.263903 0.244674 0.269559 0.250084 0.272653 0.285027 0.30398 0.22964 0.285932 0.313465 0.339638 0.340962 0.288767 0.279654 0.334098 0.191417 0.292682 0.321209 0.323623 0.21826 0.306925 0.281739 0.156824 0.307748 0.280431 0.214698 0.251721 0.356498 0.352179 0.172673 0.349373 0.352968 0.361558 0.300567 0.278045 0.30444 0.297318 0.310165 0.237732 0.251125 0.358356 0.142642 0.334533 0.249347 0.360509 0.292221 0.360263 0.210976 0.311825 0.330302 0.291789 0.289201 0.313793 0.302327 0.32703 0.252919 0.224228 0.339626 0.219558 0.296253 0.249119 0.242855 0.237096 0.236922 0.278599 0.249532 0.342489 0.359519 0.187167 0.30753 0.3739 0.270738 0.358625 0.245182 0.240315 0.268952 0.283816 0.208431 0.311491 0.294801 0.342395 0.293778 0.222283 0.283954 0.277997 0.367379 0.275439 0.276112 0.306636 0.239455 0.32961 0.261615 0.294179 0.235058 0.321259 0.299774 0.362337 0.360775 0.315729 0.206431 0.309134 0.250141 0.215856 0.281359 0.294903 0.282096 0.350079 0.069653 0.339828 0.196725 0.294571 0.30936 0.280295 0.286947 0.104992 0.269578 0.320464 0.285105 0.264115 0.231896 0.300979 0.321771 0.300899 0.288347 0.382012 0.293387 0.343614 0.289369 0.131241 0.201354 0.345385 -0.00078 0.342952 0.319711 0.225117 0.298114 0.33114 0.292429 0.2904 0.31141 0.264357 0.269819 0.264827 0.273722 0.258862 0.441896 0.284621 0.349058 0.350249 0.255424 0.284683 0.296648 0.359404 0.254408 0.296888 0.297378 0.20215 0.309002 0.253298 0.359091 0.242188 0.264315 0.308096 0.284992 0.07294 0.348489 0.243821 0.286116 0.244578 0.270672 0.226653 0.287426 0.210405 0.319221 0.246903 0.379613 0.291554 0.311444 0.325002 0.201174 0.365433 0.158248 0.309042 0.249132 0.184367 0.291351 0.319204 0.281069 0.228828 0.377707 0.310127 0.269354 0.256579 0.28595 0.2801 0.282112 0.253075 0.339099 0.346039 0.278158 0.438364 0.230964 0.159379 0.325802 0.143771 0.261995 0.337244 0.270492 0.277977 0.355962 0.379028 0.321361 0.321304 0.298468 0.262494 0.31018 0.317571 0.31743 0.299176 0.22737 0.177551 0.336749 0.280514 0.202252 0.35401 0.261842 0.331171 0.231012 0.347279 0.204604 0.319956 0.3202 0.394385 0.243243 0.284166 0.235175 0.312272 0.34459 0.315408 0.296123 0.242228 0.219202 0.267158 0.281193 0.298244 0.193979 0.375427 0.252544 0.26588 0.238385 0.291589 0.326967 0.319488 0.297156 0.282937 0.255888 0.309691 0.157182 0.362476 0.326389 0.367992 0.273281 0.313905 0.311625 0.316153 0.228216 0.326904 0.229832 0.193854 0.260013 0.322904 0.236572 0.375788 0.264683 0.289551 0.334607 0.341135 0.32623 0.205547 0.320387 0.277112 0.320558 0.244935 0.385481 0.332662 0.269058 0.316985 0.251168 0.287657 0.254978 0.116862 0.306279 0.442311 0.173671 0.260886 0.234573 0.333486 0.334924 0.361909 0.345442 0.276218 0.213197 0.243321 0.267814 0.269128 0.300885 0.342684 0.352392 0.309854 0.303761 0.319792 0.347337 0.345828 0.242363 0.30657 0.310984 0.352949 0.204763 0.323282 0.457 0.278061 0.304843 0.301948 0.065826 0.25524 0.305652 0.231737 0.452034 0.404243 0.335598 0.250436 0.309449 0.28206 0.240044 0.285716 0.225279 0.307466 0.309004 0.363164 0.297939 0.299003 0.271968 0.367708 0.257867 0.304253 0.29376 0.317879 0.291997 0.208353 0.290041 0.127643 0.21969 0.342241 0.384799 0.354736 0.194858 0.304825 0.303706 0.343584 0.369547 0.304434 0.413726 0.274517 0.259007 0.251735 0.21857 0.26166 0.235975 0.321301 0.2794 0.257327 0.306625 0.280469 0.323608 0.205092 0.397468 0.3512 0.246983 0.299266 0.22596 0.225765 0.323637 0.305028 0.344371 0.371413 0.27045 0.28073 0.285921 0.318257 0.378848 0.321804 0.257519 0.30669 0.347967 0.309297 0.328486 0.270754 0.317124 0.313877 0.343437 0.279189 0.348761 0.296603 0.329595 0.304701 0.271024 0.228349 0.302257 0.328752 0.236025 0.256284 0.243216 0.336223 0.373056 0.184801 0.232819 0.298875 0.271248 0.224524 0.329801 0.269632 0.343902 0.360462 0.351066 0.327003 0.390651 0.343521 0.31729 0.326723 0.389756 0.292872 0.250875 0.368597 0.316325 0.239098 0.286149 0.32693 0.332039 0.317395 0.254388 0.331188 0.250551 0.3278 0.305822 0.29912 0.263118 0.337704 0.304165 0.256063 0.325391 0.237902 0.257653 0.234055 0.308551 0.254192 0.237809 0.321227 0.142999 0.244137 0.314576 0.334502 0.285847 0.356333 0.179961 0.288292 0.275237 0.283543 0.279066 0.20839 0.311532 0.311842 0.286379 0.175992 0.241415 0.466064 0.322855 0.370847 0.295726 0.190451 0.17459 0.314207 0.337485 0.372317 0.232816 0.313608 0.328395 0.257448 0.358016 0.357919 0.329035 0.330506 0.308878 0.222722 0.330571 0.304504 0.357113 0.219154 0.370012 0.201279 0.339163 0.261353 0.347027 0.303045 0.347944 0.268609 0.338211 0.283501 0.373942 0.370656 0.428569 0.251905 0.347402 0.29911 0.350979 0.261287 0.272686 0.342341 0.322122 0.297719 0.291457 0.172347 0.380812 0.073102 0.255112 0.305787 0.371071 0.274671 0.079092 0.405131 0.251763 0.288674 0.284497 0.394343 0.376551 0.196202 0.214686 0.280731 0.228305 0.232897 0.288079 0.268826 0.335862 0.270489 0.319161 0.383854 0.297618 0.359975 0.321158 0.220593 0.333557 0.321844 0.290217 0.273707 0.266492 0.319063 0.346811 0.291671 0.221503 0.313775 0.338612 0.292381 0.355366 0.319876 0.383812 0.258871 0.202779 0.241203 0.287923 0.202253 0.330577 0.235902 0.299849 0.420834 0.318475 0.291765 0.348721 0.261069 0.184283 0.309528 0.098528 0.379766 0.286709 0.304364 0.390971 0.346433 0.303143 0.313388 0.309037 0.302222 0.288081 0.285691 0.361086 0.369565 0.322486 0.229466 0.2878 0.300514 0.305393 0.150756 0.321296 0.312944 0.301799 0.254133 0.272326 0.271131 0.358957 0.283038 0.260054 0.307841 0.236861 0.135106 0.303221 0.349968 0.279672 0.320254 0.218679 0.375462 0.29757 0.256886 0.1925 0.231668 0.220634 0.306402 0.33504 0.203733 0.257777 0.217464 0.16686 0.260822 0.375577 0.303675 0.316801 0.307052 0.352875 0.330808 0.301883 0.192618 0.381873 0.374069 0.365841 0.337114 0.354327 0.274295 0.29299 0.279503 0.339118 0.348654 0.2374 0.30339 0.173736 0.236991 0.285318 0.361019 0.288473 0.309497 0.289269 0.274392 0.219441 0.229052 0.259962 0.277579 0.232936 0.253729 0.222155 0.33989 0.270037 0.297987 0.414614 0.326041 0.158168 0.250509 0.268328 0.079551 0.335231 0.316968 0.39359 0.269395 0.238999 0.254987 0.28281 0.417973 0.458084 0.328101 0.349979 0.341341 0.335202 0.322604 0.359734 0.281041 0.284521 0.277166 0.259976 0.279888 0.272278 0.312078 0.323055 0.32551 0.320085 0.429036 0.287218 0.377243 0.262693 0.321063 0.175316 0.25359 0.33478 0.256664 0.295736 0.331685 0.161748 0.333196 0.246592 0.158194 0.266083 0.09904 0.282061 0.38405 0.209847 0.281617 0.206923 0.322075 0.306253 0.26453 0.325936 0.233526 0.249671 0.272573 0.25225 0.276579 0.264354 0.320184 0.285827 0.139843 0.317884 0.17497 0.261918 0.297019 0.37635 0.277017 0.277548 0.344983 0.327979 0.291471 0.292285 0.370041 0.215746 0.34124 0.290911 0.323739 0.273374 0.325878 0.306426 0.310578 0.367255 0.31056 0.170185 0.238057 0.247344 0.285505 0.15408 0.329712 0.340132 0.276026 0.4137 0.275092 0.36341 0.253369 0.331091 0.336364 0.288159 0.238481 0.302366 0.270125 0.251129 0.276979 0.236556 0.271366 0.349301 0.23606 0.25596 0.065156 0.334494 0.391953 0.252559 0.170344 0.289044 0.234151 0.360227 0.250545 0.235642 0.300876 0.258127 0.304748 0.259196 0.245066 0.265883 0.285671 0.258588 0.265231 0.218902 0.337028 0.308043 0.275396 0.258793 0.238335 0.313586 0.354255 0.290952 0.159805 0.264707 0.311639 0.306561 0.308977 0.318125 0.30716 0.319764 0.273516 0.299859 0.311219 0.308229 0.327971 0.268202 0.33698 0.293007 0.332969 0.138986 0.308777 0.218139 0.286745 0.35683 0.346157 0.181956 0.281229 0.238776 0.276036 0.244321 0.27676 0.290509 0.298985 0.302934 0.241416 0.257095 0.360351 0.333149 0.265652 0.138126 0.364059 0.341036 0.263852 0.330484 0.266817 0.213608 0.307902 0.39785 0.325822 0.34029 0.277932 0.245487 0.195768 0.331309 0.264101 0.352497 0.246325 0.379035 0.251669 0.359877 -0.123805 0.199856 0.294648 0.34679 0.269262 0.250117 0.300574 0.314113 0.317357 0.269711 0.118809 0.30875 0.297871 0.2609 0.30468 0.283159 0.30633 0.318356 0.246612 0.357015 0.205958 0.303851 0.321746 0.324283 0.28956 0.201632 0.189304 0.180343 0.393914 0.359381 0.33154 0.268035 0.230579 0.315299 0.303507 0.331611 0.285969 0.250129 0.27189 0.223314 0.213648 0.281209 0.341373 0.009119 0.136685 0.32236 0.20913 0.219028 0.36489 0.32464 0.218964 0.382951 0.297126 0.374185 0.330048 0.265339 0.42918 0.2553 0.301458 0.269601 0.346211 0.272079 0.325533 0.375123 0.216228 0.336982 0.127204 0.324689 0.21784 0.250524 0.375128 0.254383 0.309348 0.096986 0.3016 0.340567 0.29363 0.239268 0.379271 0.342362 0.326936 0.351026 0.22667 0.355388 0.302739 0.313753 0.258243 0.312329 0.267524 0.241156 0.276949 0.230185 0.266682 0.28583 0.289187 0.287479 0.309339 0.360772 0.286716 0.328317 0.303922 0.039062 0.326924 0.255302 0.242464 0.273718 0.212656 0.305076 0.239973 0.313996 0.264741 0.272697 0.221123 0.253295 0.290261 0.371182 0.213692 0.221847 0.299911 0.291515 0.313733 0.306724 0.261285 0.381839 0.325926 0.197275 0.081551 0.268192 0.329774 0.202512 0.273553 0.307568 0.509955 0.347203 0.18844 0.335322 0.316633 0.197092 0.289538 0.285141 0.276043 0.281944 0.182636 0.317028 0.151949 0.335338 0.305037 0.312123 0.342817 0.24527 0.332151 0.377299 0.26831 0.391855 0.285083 0.280415 0.300381 0.262732 0.202164 0.237655 0.29126 0.31811 0.365804 0.225396 0.347987 0.286767 0.245179 0.282592 0.24199 0.304824 0.312945 0.327881 0.211769 0.30013 0.422746 0.258509 0.314581 0.276313 0.296546 0.310328 0.312656 0.256937 0.202763 0.414456 0.309745 0.121109 0.374706 0.361753 0.297992 0.275218 0.20108 0.309378 0.374039 0.355926 0.361466 0.235022 0.335977 0.26504 0.346085 0.404142 0.247635 0.392002 0.267049 0.300884 0.204117 0.258118 0.264311 0.354199 0.253969 0.235606 0.239145 0.241038 0.294436 0.296069 0.250755 0.280382 0.298444 0.317593 0.233237 0.076348 0.280504 0.31476 0.257801 0.356332 0.286125 0.332738 0.222879 0.400843 0.353607 0.268249 0.280847 0.238663 0.381993 0.287544 0.164771 0.275841 0.242172 0.326594 0.261649 0.2602 0.262153 0.298601 0.284403 0.279332 0.337143 0.316249 0.274054 0.33338 0.165929 0.259463 0.265449 0.115063 0.266416 0.32213 0.321519 0.276434 0.26074 0.288649 0.332606 0.329272 0.36914 0.44139 0.292243 0.247603 0.207375 0.307399 0.310577 0.269704 0.268967 0.248282 0.322588 0.267121 0.167632 0.292969 0.244803 0.289093 0.309987 0.294465 0.351333 0.351103 0.276401 0.200276 0.181903 0.284451 0.365263 0.299983 0.34649 0.454955 -0.009879 0.295294 0.383904 0.340576 0.261392 0.309015 0.311978 0.279593 0.279484 0.342185 0.304235 0.317913 0.27207 0.271505 0.201358 0.267136 0.270175 0.279025 0.277962 0.239549 0.192663 0.232031 0.309226 0.320523 0.146434 0.345782 0.300622 0.305682 0.277567 0.245242 0.409989 0.244077 0.344124 0.246036 0.251045 0.282373 0.253951 0.295512 0.210105 0.295562 0.286927 0.30473 0.344081 0.180664 0.314968 0.092298 0.350438 0.337346 0.430985 0.272035 0.232063 0.207645 0.308937 0.298778 0.326155 0.371096 0.270313 0.209843 0.190642 0.375024 0.295265 0.373051 0.196505 0.299914 0.279174 0.21446 0.317474 0.017931 0.326907 0.297392 0.176086 0.292584 0.268634 0.258468 0.295969 0.292978 0.243558 0.248069 0.288168 0.330209 0.272068 0.262964 0.267205 0.283029 0.297144 0.224483 0.274167 0.264337 0.353311 0.346822 0.309764 0.244341 0.274692 0.287417 0.316497 0.005698 0.150061 0.281987 0.312183 0.294188 0.236295 0.26517 0.324888 0.297188 0.240884 0.291115 0.367944 0.289074 0.3552 0.304158 0.308498 0.314192 0.319287 0.3409 0.121714 0.252377 0.190801 0.335684 0.295605 0.131087 0.292008 0.431026 0.201493 0.309933 0.330257 0.21362 0.299137 0.351091 0.231787 0.320639 0.390553 0.333471 0.177352 0.292075 0.350245 0.306799 0.298334 0.3734 0.120698 0.313868 0.321934 0.260306 0.289754 0.089498 0.082203 0.17381 0.309061 0.288373 0.27886 0.209107 0.299168 0.303146 0.260378 0.291488 0.206744 0.191063 0.305722 0.28965 0.284368 0.260512 0.315864 0.275035 0.079746 0.284706 0.21607 0.3141 0.276224 0.292378 0.272719 0.294283 0.275625 0.265172 0.140007 0.341279 0.337977 0.310788 0.372675 0.225697 0.270554 0.265519 0.32912 0.324517 0.259755 0.287278 0.242044 0.130906 0.29859 0.262355 0.291186 0.203495 0.305513 0.273033 0.34647 0.362402 0.326729 0.361195 0.461283 0.309771 0.214746 0.269318 0.352548 0.421625 0.329744 0.324928 0.369988 0.338806 0.271027 0.330719 0.121283 0.312349 0.30423 0.320551 0.217105 0.172386 0.262822 0.356475 0.303054 0.305638 0.192949 0.317085 0.281082 0.368009 0.334809 0.156657 0.225357 0.285159 0.170183 0.223484 0.17024 0.30115 0.270598 0.259915 0.260446 0.037044 0.269395 0.264014 0.217888 0.268552 0.180968 0.20592 0.378085 0.365428 0.309943 0.302255 0.285741 0.196707 0.262079 0.375659 0.310472 0.294584 0.331756 0.347721 0.168787 0.302241 0.144101 0.353977 0.224262 0.208332 0.234448 0.306715 0.359373 0.335057 0.339338 0.296354 0.348448 0.245911 0.364459 0.202025 0.281863 0.287916 0.317793 0.175543 0.291268 0.301091 0.29509 0.334269 0.282651 0.254435 0.244616 0.24269 0.274654 0.350043 0.299406 0.316726 0.261523 0.303013 0.252475 0.168698 0.313984 0.254257 0.311722 0.194931 0.343396 0.242244 0.363169 0.386721 0.259019 0.23826 0.437972 0.314607 0.289485 0.366983 0.276371 0.230795 0.212765 0.266172 0.275221 0.295575 0.278099 0.192828 0.377648 0.311595 0.223667 0.305899 0.316665 0.080331 0.26888 0.284312 0.207957 0.301203 0.172105 0.321833 0.274187 0.381156 0.125686 0.275214 0.278839 0.308837 0.34063 0.279654 0.300527 0.250909 0.342633 0.25089 0.209228 0.241617 0.320442 0.282083 0.26968 0.310991 0.316734 0.375248 0.381809 0.258081 0.333388 0.304121 0.219894 0.275146 0.247865 0.271446 0.324543 0.293956 0.336029 0.342628 0.280311 0.288413 0.293852 0.302871 0.271199 0.383259 0.288863 0.190578 0.265748 0.298917 0.36788 0.412906 0.438218 0.266771 0.382874 0.304317 0.331355 0.273674 0.264178 0.305671 0.251214 0.247309 0.315077 0.309435 0.296774 0.240803 0.254657 0.340156 0.326955 0.23103 0.346873 0.275553 0.195314 0.16387 0.419492 0.355249 0.244213 0.333857 0.062022 0.283721 0.083817 0.350946 0.319994 0.270742 0.278815 0.178209 0.231233 0.35788 0.278837 0.299635 0.248719 0.282273 0.295223 0.307873 0.282363 0.29121 0.348873 0.347839 0.238153 0.318799 0.273854 0.277493 0.320886 0.306964 0.192286 0.28954 0.233923 0.102957 0.352459 0.251294 0.331981 0.334427 0.3003 0.315136 0.203883 0.243437 0.181479 0.317883 0.328521 0.263171 0.281137 0.262378 0.270451 0.296011 0.37658 0.301047 0.234301 0.297666 0.221574 0.180098 0.280593 0.246674 0.393282 0.324677 0.278071 0.231653 0.233304 0.211195 0.289606 0.311518 0.205776 0.299147 0.301069 0.295714 0.273193 0.346709 0.260819 0.436203 0.284868 0.263503 0.352603 0.177138 0.164265 0.315966 0.10363 0.252044 0.199997 0.382214 0.115111 0.308824 0.189524 0.208263 0.229505 0.295143 0.301882 0.327433 0.306031 0.398328 0.275123 0.403947 0.275314 0.292377 0.223485 0.26229 0.362978 0.358375 0.369834 0.277753 0.185051 0.336182 0.122065 0.319121 0.291378 0.260766 0.330219 0.295095 0.257931 0.352479 0.192244 0.321563 0.243886 0.30693 0.323569 0.411985 0.341264 0.31422 0.28506 0.281202 0.346702 0.248833 0.363595 0.270045 0.335451 0.31905 0.320693 0.381253 0.23907 0.304826 0.245095 0.242657 0.321001 0.236975 0.270321 0.331295 0.208247 0.265394 0.330296 0.294332 0.282675 0.302645 0.287612 0.305143 0.303668 0.355517 0.245521 0.269592 0.343526 0.309351 0.235579 0.218091 0.368184 0.278935 0.214056 0.25234 0.237775 0.21287 0.340561 0.297718 0.302936 0.131211 0.243671 0.359041 0.351164 0.307398 0.241693 0.312816 0.309662 0.40422 0.317476 0.357403 0.239115 0.253463 0.234282 0.191929 0.263912 0.353823 0.314395 0.297421 0.261616 0.33599 0.341698 0.363274 0.343219 0.342628 0.228026 0.272827 0.308501 0.264695 0.27949 0.271226 0.281093 0.31838 0.321674 0.286988 0.230198 0.279007 0.336397 0.316246 0.335674 0.280622 0.265651 0.338277 0.224638 0.281778 0.24354 0.332515 0.375174 0.167017 0.349436 0.328866 0.359018 0.345128 0.315752 0.188488 0.383489 0.308862 0.342064 0.109312 0.145484 0.146981 0.300331 0.210599 0.313544 0.292669 0.323668 0.329544 0.267086 0.343732 0.304108 0.326106 0.316045 0.297493 0.25233 0.316863 0.317339 0.352456 0.313661 0.139445 0.331383 0.248194 0.264899 0.240365 0.365129 0.154092 0.311768 0.225597 0.352507 0.137988 0.243488 0.26231 0.340228 0.313962 0.325478 0.244138 0.30909 0.253891 0.217942 0.262156 0.267687 0.234407 0.076969 0.38382 0.170016 0.308195 0.311894 0.271611 0.347241 0.318725 0.222189 0.369436 0.236102 0.200973 0.296678 0.301297 0.352677 0.340875 0.270938 0.380962 0.318283 0.30448 0.335544 0.314812 0.34117 0.156443 0.150535 0.416751 0.274329 0.362414 0.401852 0.268951 0.218832 0.327026 0.292348 0.30871 0.30434 0.310928 0.252925 0.259534 0.320984 0.382988 0.368588 0.362308 0.309883 0.251919 0.250982 0.300879 0.160846 0.329168 0.053274 0.333174 0.312293 0.24622 0.240172 0.289816 0.378062 0.275315 0.264574 0.31136 0.337043 0.234476 0.325895 0.230359 0.230452 0.052605 0.242641 0.242212 0.269082 0.192746 0.324397 0.285834 0.358638 0.296251 0.281808 0.295851 0.29019 0.287323 0.217308 0.25035 0.332659 0.310502 0.273147 0.292781 0.326598 0.297286 0.194416 0.312772 0.232774 0.28859 0.300834 0.199548 0.278049 0.299577 0.343625 0.117865 0.337602 0.262569 0.44457 0.272251 0.319948 0.31867 0.319454 0.217998 0.258565 0.371872 0.301196 0.313028 0.354371 0.269255 0.415741 0.314144 0.282407 0.368047 0.30117 0.38399 0.256484 0.311776 0.317512 0.240018 0.180427 0.161484 0.284303 0.318718 0.328559 0.33018 0.2492 0.338285 0.367668 0.295412 0.332695 0.116203 0.255881 0.22703 0.262424 0.389127 0.231272 0.274144 0.247751 0.213044 0.24435 0.146897 0.269714 0.337018 0.331625 0.20228 0.33266 0.385811 0.44324 0.410103 0.344547 0.287672 0.275082 0.234341 0.241315 0.22069 0.342559 0.343366 0.194853 0.303823 0.382817 0.341765 0.214226 0.268571 0.326504 0.284877 0.178944 0.32404 0.288914 0.316069 0.236631 0.298151 0.26876 0.25568 0.313096 0.292054 0.302597 0.375022 0.296595 0.295971 0.302208 0.316115 0.252141 0.347358 0.365322 0.270571 0.179711 0.177502 0.263385 0.316383 0.211866 0.343668 0.136751 0.224293 0.32054 0.210684 0.30837 0.380537 0.306818 0.337363 0.187036 0.37814 0.067053 0.12357 0.415564 0.307823 0.270072 0.38369 0.347222 0.317433 0.260393 -0.067972 0.282558 0.354484 0.315897 0.263784 0.254708 0.27808 0.337323 0.251739 0.337264 0.290271 0.293022 0.256905 0.336727 0.291652 0.194256 0.35482 0.272059 0.333153 0.350739 0.305129 0.309985 0.258514 0.318747 0.251896 0.315517 0.311794 0.266929 0.21207 0.237312 0.289797 0.208037 0.223549 0.076341 0.327845 0.363166 0.188974 0.270001 0.186649 0.22978 0.211856 0.231052 0.353857 0.315035 0.32256 0.219733 0.206368 0.206344 0.173264 0.296341 0.288095 0.361022 0.278939 0.225849 0.354492 0.305656 0.206927 0.230183 0.202485 0.318475 0.170633 0.318155 0.287505 0.308313 0.358612 0.341361 0.265822 0.22892 0.269178 0.268787 0.25709 0.33437 0.315914 0.248668 0.334219 0.300457 0.287911 0.284022 0.342211 0.169663 0.233766 0.291024 0.306003 0.387411 0.294387 0.266379 0.312817 0.292784 0.213672 0.25763 0.265922 0.325848 0.315458 0.252394 0.311503 0.034332 0.16688 0.24982 0.265025 0.291208 0.270599 0.280678 0.232297 0.303742 0.287066 0.336879 0.294782 0.24733 0.254491 0.26115 0.452869 0.108033 0.262571 0.291204 0.256994 0.403562 0.140507 0.311363 0.235932 0.262682 0.380369 0.327396 0.166102 0.414224 0.23983 0.263042 0.224533 0.263245 0.302867 0.341442 0.294071 0.323836 0.348541 0.311656 0.199801 0.357166 0.274225 0.345537 0.294237 0.285018 0.202497 0.30656 0.219869 0.307372 0.312475 0.31268 0.250347 0.308137 0.262082 0.355927 0.251112 0.265787 0.355514 0.238976 0.257389 0.412281 0.274344 0.293786 0.329535 0.30156 0.238756 0.262316 0.303224 0.294647 0.336385 0.268868 0.307491 0.295279 0.302908 0.312113 0.31377 0.290481 0.098716 0.38256 0.280047 0.264896 0.231372 0.375022 0.271604 0.32621 0.279544 0.108779 0.324837 0.248061 0.278801 0.322764 0.352591 0.308173 0.241146 0.267898 0.202541 0.314652 0.233661 0.226968 0.280877 0.256825 0.241407 0.303466 0.21991 0.101482 0.369271 0.253708 0.331877 0.336822 0.33184 0.387858 0.315578 0.357086 0.210725 0.268427 0.313669 0.173672 0.325755 0.289151 0.302717 0.25113 0.226689 0.330057 0.324638 0.257294 0.32439 0.261154 0.275525 0.327387 0.201528 0.194498 0.275839 0.290764 0.240575 0.220991 0.324673 0.337621 0.283254 0.210273 0.278768 0.343735 0.203538 0.340177 0.276 0.329933 0.26096 0.312955 0.307745 0.341718 0.263467 0.239388 0.324007 0.289912 0.277707 0.309696 0.393845 0.210178 0.263324 0.26965 0.097856 0.264562 0.299455 0.294572 0.262463 0.261275 0.254113 0.27272 0.272609 0.237311 0.287405 0.268986 0.304525 0.211253 0.225809 0.41177 0.300196 0.149492 0.252284 0.311598 0.24058 0.39876 0.204615 0.164521 0.322962 0.376257 0.312516 0.295408 0.35023 0.193801 0.337366 0.252093 0.297039 0.320429 0.198 0.004551 0.218536 0.452555 0.400848 0.291026 0.148934 0.158435 0.315701 0.261723 0.138803 0.470376 0.257888 0.34392 0.326588 0.249272 0.334883 0.339346 0.228604 0.320655 0.276946 0.244463 0.242581 0.232536 0.323166 0.27585 0.221446 0.244599 0.383578 0.280272 0.111801 0.127774 0.254124 0.284476 0.305497 0.275533 0.293519 0.32866 0.340541 0.278447 0.211636 0.239016 0.337951 0.295908 0.304295 0.339333 0.319499 0.280706 0.307024 0.281067 0.311685 0.325908 0.297633 0.33788 0.201622 0.244394 0.303176 0.264479 0.232542 0.319227 0.364796 0.204217 0.292138 0.274539 0.354434 0.246343 0.257277 0.237482 0.181522 0.166678 0.290007 0.364247 0.311947 0.2265 0.258968 0.273751 0.23587 0.298294 0.235691 0.273878 0.369692 0.073069 0.381378 0.22152 0.260242 0.33935 0.092409 0.325313 0.330835 0.276621 0.361317 0.267947 0.194851 0.313978 0.323821 0.35229 0.224272 0.305305 0.215192 0.273898 0.262956 0.333263 0.309242 0.250786 0.315041 0.378272 0.372983 0.337442 0.20007 0.325782 0.298212 0.333755 0.439501 0.371939 0.288617 0.264502 0.186606 0.322578 0.308727 0.252809 0.343079 0.359382 0.204451 0.366286 0.301553 0.309173 0.27177 0.2733 0.19852 0.281795 0.30794 0.342307 0.303619 0.364576 0.32955 0.307702 0.268969 0.317854 0.181793 0.205351 0.238968 0.278232 0.30888 0.319645 0.2235 0.257667 0.285462 0.34257 0.251572 0.289066 0.265287 0.252887 0.352023 0.127326 0.359277 0.173942 0.277209 0.307177 0.254162 0.341858 0.270647 0.277194 0.458693 0.217093 0.304503 0.306351 0.250465 0.289413 0.254423 0.279287 0.401463 0.328061 0.364994 0.305463 0.268584 0.320855 0.322973 0.329055 0.240735 0.279534 0.195155 0.333137 0.222972 0.255765 0.312275 0.267896 0.314268 0.293862 0.260555 0.22757 0.333429 0.296199 0.269663 0.291094 0.317512 0.317039 0.369463 0.293692 0.382621 0.226523 0.271424 0.267162 0.382799 0.18356 0.045583 0.262196 0.339121 0.305458 0.33677 0.336417 0.260236 0.259222 0.34593 0.224658 0.330632 0.30765 0.319739 0.297766 0.275959 0.257115 0.333909 0.290612 0.300604 0.289809 0.275317 0.32359 0.356708 0.306831 0.288015 0.384882 0.227828 0.18774 0.289783 0.269523 0.2732 0.243405 0.334839 0.241496 0.214738 0.293775 0.302517 0.247109 0.162083 0.322979 0.23198 0.297119 0.39646 0.252908 0.371153 0.413647 0.3396 0.334875 0.298655 0.33376 0.274341 0.324338 0.251356 0.268835 0.295507 0.163688 0.331441 0.295377 0.300399 0.265291 0.410165 0.372995 0.33254 0.31845 0.36765 0.281428 0.263752 0.356964 0.291959 0.328611 0.238135 0.237742 0.349971 0.316908 0.210146 0.292437 0.249949 0.336336 0.308617 0.320934 0.343606 0.270652 0.34613 0.25038 0.381561 0.309318 0.26898 0.284219 0.296309 0.358515 0.276846 0.283158 0.265116 0.22109 0.262344 0.221017 0.105005 0.223805 0.255695 0.27552 0.350831 0.324316 0.220036 0.238671 0.330145 0.2761 0.317702 0.268759 0.289799 0.135431 0.326651 0.34387 0.323037 0.286625 0.293495 0.340319 0.279549 0.262311 0.111518 0.29583 0.323985 0.242783 0.155855 0.328274 0.251872 0.186161 0.366378 0.403942 0.325973 0.404486 0.286307 0.325992 0.224717 0.271852 0.283034 0.301825 0.268548 0.234161 0.298117 0.284987 0.220304 0.160532 0.232287 0.316069 0.279697 0.201068 0.322071 0.244493 0.289755 0.279273 0.287382 0.257812 0.287754 0.289341 0.192075 0.302975 0.328263 0.19587 0.175206 0.217006 0.276368 0.20591 0.261833 0.29889 0.234842 0.296541 0.279429 0.310231 0.277738 0.297001 0.313039 0.257555 0.22635 0.414001 0.23827 0.26942 0.325132 0.32396 0.31067 0.259109 0.294228 0.205795 0.280428 0.124496 0.267001 0.298009 0.241135 0.257652 0.21084 0.260402 0.312671 0.303045 0.265885 0.25119 0.322007 0.287132 0.234085 0.381758 0.329869 0.334179 0.253104 0.276098 0.303672 0.261663 0.441188 0.315688 0.36043 0.209664 0.283402 0.269092 0.324792 0.322035 0.299862 0.324788 0.284777 0.341217 0.252349 0.294225 0.26987 0.256677 0.288804 0.324514 0.304113 0.402565 0.336161 0.382718 0.206528 0.327264 0.272117 0.249715 0.254342 0.284387 0.352899 0.297729 0.305465 0.297745 0.175738 0.24914 0.260678 0.348799 0.112318 0.361471 0.126175 0.35145 0.201827 0.32179 0.306946 0.148366 0.303021 0.308919 0.325543 0.30782 0.273093 0.326607 0.477994 0.395075 0.262549 0.246621 0.268749 0.329863 0.364584 0.290605 0.217895 0.21791 0.215011 0.369205 0.207528 0.329901 0.251921 0.30864 0.308228 0.264875 0.320027 0.315635 0.306562 0.390365 0.265798 0.041156 0.359063 0.295043 0.26882 0.273627 0.209001 0.3246 0.165687 0.355796 0.213638 0.354093 0.30102 0.362055 0.139884 0.300429 0.378569 0.27976 0.321367 0.227468 0.242197 0.213676 0.289692 0.193606 0.193246 0.252781 0.330909 0.303775 0.263887 0.344107 0.329486 0.342742 0.296939 0.33817 0.217274 0.2839 0.311246 0.286623 0.170571 0.252266 0.342625 0.330228 0.274805 0.252395 0.221923 0.318465 0.275522 0.327478 0.299615 0.341273 0.323737 0.206242 0.158399 0.337687 0.368183 0.300705 0.276739 0.216897 0.351723 0.245392 0.276315 0.304509 0.217091 0.319864 0.34297 0.319024 0.124497 0.190559 0.210036 0.312589 0.309098 0.335331 0.266944 0.325399 0.197036 0.281955 0.264977 0.306369 0.270867 0.25165 0.293256 0.289101 0.367681 0.247715 0.30647 0.362823 0.24945 0.288873 0.377657 0.309693 0.302392 0.302339 0.301719 0.259642 0.288877 0.315779 0.320139 0.281801 0.239005 0.285173 0.302748 0.300589 0.260286 0.29771 0.240461 0.207986 0.29619 0.302169 0.284018 0.268349 0.173417 0.308067 0.19375 0.172848 0.326174 0.341344 0.388762 0.263983 0.250495 0.372518 0.316274 0.308933 0.284238 0.377184 0.248647 0.255772 0.353109 0.218008 0.34807 0.337949 0.275371 0.331678 0.331043 0.315372 0.331376 0.272255 0.224106 0.352193 0.209073 0.322758 0.343168 0.281632 0.32023 0.387856 0.05191 0.250514 0.275545 0.289582 0.290531 0.210301 0.267622 0.364384 0.304903 0.263001 0.394311 0.303112 0.320471 0.3048 0.286022 0.028837 0.3901 0.322834 0.326971 0.274572 0.289101 0.26012 0.371479 0.351763 0.324376 0.25741 0.314181 0.247575 0.356796 0.336134 0.24516 0.286684 0.247062 0.245273 0.246624 0.305159 0.36598 0.268796 0.303155 0.275178 0.332511 0.306798 0.231126 0.225322 0.265604 0.135422 0.312386 0.231861 0.260937 0.288535 0.257502 0.304116 0.260979 0.297874 0.303347 0.162418 0.261109 0.205688 0.322154 0.359174 0.329382 0.329788 0.142257 0.212672 0.290349 0.303121 0.254684 0.200802 0.283996 0.214801 0.31173 0.275361 0.177592 0.283699 0.238271 0.273069 0.241287 0.245343 0.296585 0.289011 0.311161 0.264363 0.287342 0.188871 0.34331 0.308376 0.400042 0.35601 0.390763 0.265672 0.283429 0.197651 0.313251 0.248503 0.260996 0.240408 0.269602 0.326134 0.235572 0.347646 0.266723 0.316113 0.301269 0.298296 0.411047 0.338956 0.27797 0.314616 0.311486 0.262053 0.366501 0.281051 0.338643 0.30033 0.440275 0.267978 0.302868 0.383921 0.247251 0.31119 0.29864 0.270255 0.360176 0.314986 0.174013 0.333855 0.194994 0.231202 0.461085 0.045281 0.405896 0.309562 0.288449 0.251264 0.306858 0.239051 0.229162 0.231566 0.215673 0.361728 0.295484 0.251351 0.26875 0.319158 0.342217 0.295945 0.317139 0.272786 0.251372 0.331189 0.156852 0.234354 0.346979 0.272502 0.344535 0.351106 0.307142 0.382713 0.289596 0.305506 0.348118 0.195762 0.307859 0.291968 0.349627 0.153291 0.273802 0.299536 0.255535 0.255883 0.303186 0.344089 0.313036 0.31022 0.261384 0.277489 0.225393 0.313483 0.279725 0.279803 0.343398 0.257936 0.198352 0.328961 0.1948 0.189524 0.304568 0.315819 0.273924 0.161669 0.238597 0.230528 0.321657 0.292201 0.284007 0.278379 0.22944 0.291168 0.353261 0.313522 0.268966 0.298382 0.214539 0.327139 0.264458 0.4125 0.163931 0.177524 0.281138 0.329249 0.323923 0.3526 0.240285 0.317589 0.296649 0.311109 0.288871 0.341771 0.350369 0.288976 0.295883 0.327155 0.33471 0.285381 0.254258 0.302152 0.283071 0.257484 0.164227 0.327654 0.271833 0.253085 0.214694 0.205908 0.064635 0.27009 0.32664 0.316629 0.392628 0.332501 0.327084 0.338373 0.349885 0.283044 0.292658 0.274377 0.328454 0.265978 0.287009 -0.046259 0.234012 0.283411 0.276423 0.350158 0.326679 0.255516 0.297336 0.26761 0.316271 0.264873 0.32032 0.260274 0.237339 0.205567 0.28632 0.269141 0.300585 0.301452 0.319576 0.306163 0.246685 0.31255 0.255778 0.234536 0.356739 0.297017 0.303807 0.265783 0.500088 0.236164 0.222315 0.235792 0.194967 0.28365 0.3618 0.253546 0.234574 0.259881 0.267344 0.285407 0.224173 0.31847 0.331752 0.289061 0.263646 0.331031 0.328688 0.347313 0.22787 0.252849 0.247429 0.288803 0.316622 0.272869 0.315794 0.168017 0.335759 0.274006 0.261333 0.213114 0.14822 0.349232 0.340014 0.314225 0.29172 0.300085 0.333613 0.314471 0.278105 0.250729 0.2978 0.206006 0.398041 0.227279 0.263031 0.275401 0.297708 0.317769 0.253143 0.297123 0.316121 0.260147 0.197738 0.415922 0.356883 0.266747 0.246825 0.309759 0.263299 0.344411 0.242514 0.385619 0.281219 0.361151 0.415148 0.313074 0.316109 0.204305 0.297969 0.295264 0.057229 0.243222 0.31295 0.311213 0.273077 0.2531 0.242279 0.295685 0.331693 0.280066 0.357744 0.241554 0.083407 0.271246 0.233069 0.322366 0.389929 0.34959 0.443639 0.227808 0.242118 0.352175 0.339318 0.397222 0.308582 0.188828 0.278451 0.280009 0.159175 0.286364 0.282597 0.228692 0.296243 0.384179 0.290601 0.334302 0.26244 0.304765 0.298908 0.275525 0.180724 0.246622 0.19879 0.326873 0.279096 0.225922 0.263502 0.296834 0.285189 0.29846 0.41138 0.334132 0.303393 0.261995 0.232301 0.081086 0.305084 0.334991 0.276604 0.268913 0.340531 0.254426 0.346359 0.243058 0.305267 0.365261 0.294369 0.324366 0.240205 0.237366 0.175062 0.302455 0.284079 0.348622 0.363699 0.243566 0.27304 0.289619 0.19545 0.325533 0.320697 0.291743 0.345202 0.244007 0.278025 0.264047 0.275228 0.309339 0.351766 0.360838 0.312016 0.289294 0.348166 0.300692 0.23414 0.311161 0.148444 0.293493 0.298631 0.320351 0.258086 0.270228 0.348186 0.371707 0.396266 0.321308 0.250749 0.346932 0.283949 0.303393 0.367405 0.265818 0.241941 0.314447 0.32543 0.203683 0.30873 0.220873 -0.036953 0.262053 0.259935 0.202279 0.158871 0.31186 0.335134 -0.065117 0.334211 0.236459 0.175253 0.247296 0.25784 0.368802 0.231598 0.250289 0.284162 0.352938 0.304124 0.24863 0.310363 0.377741 0.293074 0.326713 0.334417 0.382356 0.251065 0.337042 0.261355 0.111129 0.294038 0.241803 0.279432 0.289913 0.253142 0.291159 0.326796 0.135672 0.301453 0.276986 0.299801 0.203807 0.160452 0.341663 0.285358 0.281048 0.388889 0.309895 0.339343 0.198482 0.313082 0.291701 0.354164 0.287617 0.338768 0.359314 0.364128 0.360108 0.380883 0.248436 0.316999 0.351801 0.325338 0.330233 0.343336 0.33903 0.230394 0.402685 0.238489 0.255244 0.250379 0.333114 0.215936 0.274789 0.387917 0.265461 0.259925 0.306631 0.174529 0.297186 0.31655 0.289672 0.344265 0.235576 0.299197 0.208741 0.360551 0.349678 0.222004 0.326568 0.329433 0.379782 0.336298 0.285579 0.308371 0.293149 0.296459 0.251926 0.357765 0.193293 0.416782 0.180826 0.249038 0.311243 0.279609 0.210817 0.392959 0.314126 0.287217 0.334024 0.247583 0.120118 0.253983 0.272878 0.288889 0.298749 0.263469 0.200161 0.251938 0.220612 0.210327 0.323973 0.297963 0.263016 0.23925 0.330864 0.290425 0.325059 0.354587 0.320403 0.19546 0.202021 0.275621 0.284862 0.343424 0.2312 0.275595 0.123779 0.249074 0.228079 0.299069 0.272339 0.32415 0.324149 0.224601 0.409946 0.380219 0.335387 0.365645 0.27864 0.308184 0.418996 0.304972 0.291959 0.401736 0.283578 -0.027888 0.246499 0.311825 0.291153 0.33143 0.251802 0.187104 0.245969 0.341322 0.322983 0.478512 0.308737 0.289604 0.210114 0.273804 0.375351 0.306252 0.219838 0.358475 0.254833 0.383524 0.350293 0.210409 0.284945 0.289113 0.210793 0.378827 0.334284 0.272472 0.260715 0.313073 0.280626 0.22261 0.315499 0.338412 0.314506 0.342207 0.283859 0.277084 0.310801 0.272114 0.218241 0.279345 0.187474 0.279745 0.226273 0.362084 0.308814 0.273753 0.286878 0.216582 0.344357 0.274729 0.302032 0.270732 0.31563 0.147151 0.262928 0.216395 0.231898 0.324788 0.237783 0.182297 0.215482 0.289487 0.208186 0.321359 0.344835 0.205795 0.330535 0.250938 0.319931 0.399595 0.287899 0.386334 0.334968 0.314586 0.305611 0.244903 0.289253 0.370888 0.322066 0.26108 0.116115 0.352241 0.231811 0.305815 0.239157 0.347012 0.332243 0.315487 0.293912 0.314189 0.199827 0.18292 0.311722 0.279547 0.094301 0.280024 0.260363 0.356105 0.307942 0.331762 0.276431 0.254471 0.166604 0.223291 0.259066 0.306189 0.350112 0.282451 0.332929 0.330672 0.301418 0.261215 0.298445 0.273367 0.330964 0.301313 0.32589 0.319026 0.188665 0.251854 0.294383 0.361923 0.357829 0.257953 0.275493 0.252687 0.269664 0.252995 0.34836 0.312243 0.289439 0.248567 0.267553 0.268372 0.255871 0.340261 0.18257 0.236439 0.341412 0.270751 0.215034 0.240007 0.283355 0.333454 0.338996 0.279632 0.342355 0.319534 0.29603 0.319143 0.275766 0.330388 0.305281 0.183763 0.25095 0.248239 0.394461 0.343372 0.407648 0.22625 0.275965 0.321739 0.23118 0.282095 0.304238 0.292796 0.299863 0.185751 0.268353 0.292306 0.398841 0.337751 0.259842 0.286143 0.268443 0.353353 0.183015 0.339222 0.339604 0.367938 0.36081 0.282329 0.43751 0.183369 0.199236 0.230503 0.278679 0.338238 0.252767 0.324831 0.185772 0.411534 0.31949 0.340791 0.326839 0.279305 0.262937 0.253125 0.305081 0.313482 0.280924 0.341767 0.226871 0.252886 0.266725 0.281153 0.230702 0.211657 0.315601 0.209767 0.194775 0.305663 0.326219 0.342131 0.220738 0.283903 0.285012 0.209957 0.390047 0.388102 0.292477 0.324911 0.308652 0.198254 0.294973 0.224502 0.245471 0.305163 0.368421 0.287366 0.299436 0.314298 0.254797 0.296494 0.231729 0.343441 0.379275 0.161067 0.280894 0.319219 0.307658 0.225725 0.271661 0.259275 0.322541 0.264871 0.374569 0.322775 0.365504 0.219605 0.301619 0.261166 0.275197 0.33534 0.284755 0.290605 0.356524 0.252784 0.368156 0.302597 0.402356 0.272866 0.252848 0.29504 0.341169 0.361373 0.397878 0.250668 0.369531 0.253161 0.192229 0.237638 0.378148 0.210357 0.224608 0.278361 0.308442 0.329279 0.307855 0.233131 0.264461 0.317685 0.309109 0.31278 0.289068 0.33076 0.303359 0.302614 0.218074 0.323389 0.315447 0.257329 0.265526 0.192823 0.424899 0.182874 0.391823 0.375927 0.282691 0.288218 0.249718 0.323865 0.277397 0.319835 0.226411 0.265851 0.316883 0.310675 0.316294 0.238248 0.309928 0.249977 0.372383 0.217915 0.365742 0.154481 0.351153 0.343771 0.425662 0.316316 0.064222 0.197213 0.150938 0.215855 0.332505 0.356513 0.302261 0.3023 0.260128 0.23716 0.23007 0.386014 0.274398 0.288446 0.313462 0.279653 0.311552 0.148593 0.212133 0.298169 0.348792 0.316126 0.279363 0.243404 0.258234 0.372967 0.350692 0.318701 0.213418 0.279031 0.103798 0.28879 0.32864 0.327645 0.200977 0.223292 0.267299 0.235651 0.342501 0.312618 0.217723 0.298717 0.274316 0.307515 0.283003 0.319767 0.440586 0.304733 0.151557 0.17251 0.270444 0.285514 0.360174 0.281237 0.260999 0.336843 0.18583 0.305039 0.299112 0.303594 0.243714 0.338939 0.285321 0.318314 0.302886 0.30669 0.169643 0.316024 0.339898 0.292611 0.134079 0.196372 0.34773 0.232042 0.264954 0.222546 0.21892 0.321264 0.410245 0.256669 0.305036 0.226913 0.378427 0.321442 0.23348 0.301361 0.341495 0.325984 0.373447 0.264513 0.3117 0.275228 0.245102 0.253332 0.28621 0.30852 0.309977 0.334441 0.311414 0.318809 0.382441 0.248527 0.076322 0.261847 0.283952 0.30691 0.271796 0.198479 0.271207 0.29162 0.267524 0.343378 0.043631 0.326847 0.227337 0.249559 0.088765 0.273798 0.393085 0.235784 0.295664 0.345521 0.283984 0.273982 0.33457 0.04237 0.246346 0.111144 0.354415 0.235684 0.067433 0.325064 0.290773 0.286661 0.332502 0.23795 0.277242 0.188766 0.260895 0.309552 0.26621 0.360357 0.237388 0.267971 0.372418 0.031444 0.395154 0.427745 0.290231 0.304149 0.277333 0.339708 0.308869 0.192016 0.348097 0.250898 0.386863 0.264865 0.274678 0.265278 0.249164 0.289168 0.255859 0.256833 0.379016 0.294932 0.19429 0.26229 0.342017 0.338612 0.260965 0.188569 0.266573 0.302685 0.292607 0.13508 0.322124 0.281085 0.305799 0.352658 0.271383 0.237167 0.30362 0.417706 0.292259 0.267287 0.298214 0.257544 0.305915 0.342361 0.280423 0.24835 0.324121 0.291786 0.349556 0.083873 0.312974 0.212937 0.29678 0.313323 0.318353 0.274505 0.204847 0.275486 0.292264 0.264843 0.331136 0.232711 0.313042 0.283566 0.320878 0.253795 0.322789 0.306298 0.2436 0.196008 0.252145 0.322161 0.253391 0.356175 0.242841 0.246996 0.266466 0.434169 0.312422 0.34664 0.301526 0.354858 0.234317 0.279814 0.291286 0.411671 0.278646 0.271611 0.245095 0.288117 0.407341 0.263025 0.284351 0.249793 0.267218 0.322126 0.266956 0.262463 0.301096 0.273414 0.34973 0.24191 0.294372 0.37409 0.227124 0.333374 0.303341 0.293978 0.203126 0.247444 0.205688 0.360453 0.271349 0.287088 0.303125 0.240962 0.229877 0.337145 0.313519 0.237523 0.289846 0.320493 0.343254 0.298168 0.293931 0.310886 0.151486 0.188586 0.328603 0.327683 0.282014 0.292183 0.262149 0.287644 0.281871 0.27379 0.236759 0.270465 0.337064 0.256419 0.315131 0.31907 0.346132 0.30287 0.112122 0.320892 0.265116 0.286978 0.287584 0.263182 0.257721 0.275871 0.268689 0.327268 0.326905 0.333985 0.353448 0.262209 0.296226 0.356049 0.242825 0.352845 0.375325 0.031275 0.317372 0.310068 0.289509 0.368637 0.271283 0.321195 0.36289 0.300524 0.240804 0.092962 0.225777 0.310847 0.28211 0.353533 0.247478 0.269774 0.259822 0.20845 0.311244 0.278329 0.213595 0.291097 0.259461 0.217411 0.311606 0.278511 0.288404 0.277509 0.338205 0.370016 0.354705 0.319497 0.325475 0.329939 0.381924 0.428718 0.186864 0.337401 0.350323 0.238285 0.390795 0.351917 0.294766 0.190499 0.259823 0.278461 0.09114 0.276482 0.309916 0.214918 0.27844 0.343543 0.319555 0.240629 0.240414 0.394524 0.293173 0.175213 0.268104 0.282356 0.335235 0.323242 0.154736 0.438233 0.333001 0.360577 0.349394 0.291505 0.298774 0.301939 0.141883 0.366233 0.173217 0.296275 0.157068 0.27112 0.281945 0.346701 -0.077297 0.359722 0.231982 0.266838 0.234011 0.268653 0.349809 0.270065 0.334016 0.308525 0.201684 0.326615 0.283336 0.329346 0.372927 0.284523 0.265837 0.229443 0.318345 0.293612 0.3219 0.340368 0.250474 0.288722 0.318589 0.267011 0.315918 0.377611 0.276072 0.308429 0.324236 0.354576 0.180917 0.298826 0.30291 0.347507 0.357972 0.289545 0.319334 0.329384 0.344179 0.157053 0.290619 0.295168 0.269757 0.310445 0.295258 0.307937 0.232606 0.298955 0.192413 0.225778 0.241556 0.273318 0.333553 0.256773 0.307491 0.25568 0.215605 0.343322 0.215717 0.356452 0.354215 0.332444 0.373654 0.221416 0.309105 0.23773 0.334438 0.3088 0.263851 0.274876 0.116816 0.249169 0.270343 0.172128 0.281524 0.297733 0.4973 0.286634 0.291597 0.400907 0.225555 0.268755 0.293351 0.290382 0.278356 0.311605 0.279179 0.355224 0.333418 0.264998 0.254412 0.293576 0.304937 0.334013 0.296495 0.2454 0.298425 0.31091 0.37553 0.356003 0.206851 0.198006 0.34024 0.364185 0.249664 0.333647 0.330164 0.311958 0.39172 0.238925 0.174649 0.18509 0.087011 0.310381 0.254592 0.282831 0.322971 0.275075 0.278772 0.323198 0.29594 0.265047 0.202317 0.365584 0.364915 0.303456 0.237704 0.26144 0.209296 0.299144
Exponential_normal 0.200173 0.21216 0.240657 0.224744 0.232801 0.170833 0.359593 0.16823 0.316054 0.163514 0.185518 0.162672 0.149088 0.172953 0.261438 0.160444 0.28751 0.185454 0.235242 0.247567 0.275238 0.332979 0.318688 0.201783 0.261726 0.254006 0.224616 0.176666 0.213961 0.329717 0.352177 0.195221 0.230814 0.269708 0.195116 0.328739 0.155615 0.315688 0.190253 0.184642 0.269052 0.214284 0.291117 0.219969 0.186079 0.152547 0.230101 0.184715 0.270791 0.246879 0.260523 0.257178 0.19714 0.263383 0.187829 0.17485 0.322841 0.222894 0.418057 0.212804 0.160645 0.252432 0.232835 0.237105 0.272469 0.214812 0.185953 0.219064 0.214275 0.152992 0.248199 0.237586 0.261345 0.292961 0.209329 0.25772 0.270388 0.162262 0.216891 0.326754 0.225521 0.188138 0.244202 0.210606 0.17662 0.240544 0.213487 0.157409 0.224128 0.239693 0.15838 0.232256 0.313762 0.235167 0.256262 0.233781 0.215033 0.209023 0.151556 0.326665 0.165698 0.165549 0.202704 0.244977 0.206313 0.161891 0.170003 0.247776 0.215834 0.206943 0.293861 0.189359 0.28732 0.200247 0.219586 0.223725 0.163136 0.252953 0.27171 0.198796 0.220788 0.200261 0.171423 0.161498 0.173567 0.208437 0.223558 0.182991 0.161859 0.177864 0.200898 0.183482 0.197067 0.199196 0.23396 0.177338 0.271736 0.188296 0.234333 0.167979 0.253351 0.129947 0.384427 0.188576 0.168525 0.268799 0.213996 0.187641 0.150695 0.28541 0.14473 0.250261 0.13671 0.181147 0.208771 0.188801 0.179722 0.192041 0.217457 0.234596 0.176805 0.161223 0.198146 0.183253 0.184991 0.201152 0.196573 0.125849 0.213982 0.285842 0.172971 0.208592 0.252941 0.199629 0.181058 0.133299 0.205636 0.225838 0.195902 0.152691 0.261762 0.1665 0.221894 0.24143 0.165065 0.143844 0.225022 0.238681 0.262075 0.254863 0.12525 0.241717 0.193458 0.212038 0.235574 0.17157 0.166742 0.176507 0.187163 0.30509 0.310268 0.186681 0.179948 0.190189 0.234777 0.20024 0.181793 0.20373 0.219111 0.322033 0.216406 0.198666 0.19613 0.260977 0.180558 0.177308 0.201174 0.266845 0.262448 0.204402 0.167605 0.146327 0.272251 0.18995 0.178861 0.251887 0.219046 0.177863 0.295373 0.19699 0.29523 0.311504 0.170148 0.150483 0.320074 0.219575 0.320733 0.204057 0.282152 0.181568 0.439984 0.210212 0.214895 0.214833 0.208826 0.300889 0.195243 0.232806 0.258582 0.227263 0.258442 0.176246 0.155405 0.223077 0.145817 0.155701 0.244727 0.208316 0.186791 0.201458 0.195419 0.285701 0.160844 0.21653 0.390705 0.195132 0.203295 0.23244 0.169522 0.181421 0.201634 0.174224 0.146961 0.13413 0.241311 0.279477 0.208782 0.285722 0.235652 0.15982 0.195624 0.141424 0.201849 0.223883 0.219776 0.175976 0.169581 0.215241 0.234833 0.252198 0.320364 0.221533 0.267456 0.139774 0.312407 0.187529 0.28745 0.202861 0.177731 0.294961 0.344282 0.230781 0.179414 0.152343 0.185279 0.178052 0.168122 0.22286 0.289798 0.197746 0.301583 0.167741 0.183405 0.231829 0.275458 0.279757 0.152175 0.18903 0.307469 0.201174 0.176097 0.29027 0.287245 0.197053 0.266656 0.20843 0.210989 0.184621 0.157975 0.233372 0.327889 0.270561 0.183107 0.249635 0.288134 0.210176 0.244662 0.172069 0.171534 0.167997 0.26633 0.208142 0.181211 0.275191 0.220295 0.192279 0.211241 0.160784 0.227467 0.171458 0.341327 0.210653 0.187063 0.193666 0.209291 0.169729 0.194839 0.166106 0.26725 0.205731 0.222087 0.199318 0.297074 0.26538 0.202183 0.176171 0.168002 0.262331 0.180569 0.297022 0.240079 0.206211 0.167295 0.147273 0.35743 0.189353 0.262407 0.193436 0.228451 0.172424 0.217019 0.218828 0.215594 0.203071 0.189742 0.180869 0.169098 0.187502 0.202682 0.146982 0.226725 0.220335 0.208908 0.216698 0.183645 0.186679 0.427477 0.259659 0.197872 0.168291 0.26054 0.240209 0.317122 0.229486 0.184247 0.206484 0.367607 0.232431 0.283101 0.225498 0.157244 0.395597 0.228007 0.231995 0.20283 0.249465 0.285046 0.170923 0.23148 0.221436 0.139667 0.352076 0.192406 0.158568 0.282826 0.213845 0.224222 0.262759 0.240257 0.225776 0.204712 0.299048 0.174627 0.214467 0.185974 0.248357 0.179187 0.173253 0.154904 0.204721 0.395494 0.202005 0.264022 0.198781 0.140635 0.234181 0.125429 0.256222 0.170544 0.176058 0.232913 0.17981 0.168107 0.213589 0.284865 0.229723 0.30634 0.51013 0.187527 0.227323 0.209966 0.168607 0.2101 0.232746 0.258785 0.212881 0.209587 0.216785 0.216606 0.197103 0.204809 0.142996 0.200179 0.156934 0.173578 0.167977 0.226823 0.304992 0.203202 0.160645 0.170646 0.233767 0.196748 0.157624 0.218322 0.179328 0.17584 0.237097 0.155734 0.282303 0.14008 0.306066 0.221746 0.238415 0.253655 0.185181 0.250171 0.149592 0.202103 0.229428 0.130573 0.167353 0.232461 0.217098 0.231195 0.219781 0.182608 0.299151 0.282091 0.153191 0.246882 0.225269 0.332605 0.230214 0.239651 0.212631 0.399378 0.266558 0.253155 0.166897 0.219907 0.255787 0.28354 0.189605 0.18698 0.219733 0.234417 0.195294 0.18396 0.136705 0.130365 0.177822 0.23984 0.185624 0.21785 0.286533 0.25352 0.33384 0.214717 0.274956 0.185146 0.181907 0.203406 0.268952 0.208431 0.315579 0.210486 0.193736 0.140153 0.201457 0.183703 0.184144 0.175016 0.201064 0.217596 0.183774 0.247323 0.246205 0.178432 0.22258 0.172732 0.283737 0.322457 0.177446 0.218573 0.126839 0.19674 0.197136 0.243039 0.204931 0.222264 0.186391 0.243584 0.248005 0.186405 0.21256 0.220281 0.225539 0.152236 0.237319 0.145485 0.271756 0.238996 0.225119 0.136059 0.163832 0.349786 0.181607 0.207059 0.182804 0.21106 0.179666 0.227963 0.202376 0.238165 0.191662 0.430948 0.232506 0.219789 0.17589 0.213959 0.199425 0.229911 0.181524 0.180519 0.129011 0.206865 0.168157 0.244188 0.196136 0.186225 0.170059 0.182261 0.172732 0.198232 0.279156 0.230421 0.187142 0.190584 0.208093 0.198032 0.153697 0.221375 0.290365 0.268646 0.404244 0.197941 0.155213 0.300515 0.198593 0.219246 0.139033 0.336418 0.192211 0.198324 0.20921 0.18164 0.21003 0.125691 0.187059 0.219312 0.212811 0.246793 0.179236 0.230392 0.248723 0.119025 0.204241 0.332228 0.256428 0.195662 0.171037 0.189968 0.368842 0.213623 0.25131 0.194896 0.134987 0.209229 0.226004 0.233984 0.265845 0.2035 0.18774 0.225502 0.209224 0.222939 0.187857 0.219827 0.127392 0.241898 0.123146 0.245895 0.272339 0.296851 0.181758 0.163102 0.284284 0.251228 0.181597 0.201534 0.216048 0.188258 0.134689 0.1587 0.171591 0.243533 0.143253 0.349704 0.22049 0.231079 0.140626 0.188795 0.253806 0.213824 0.335587 0.201535 0.22622 0.245846 0.378675 0.200092 0.252803 0.262484 0.287196 0.261609 0.350404 0.291926 0.225213 0.193426 0.185443 0.255973 0.129106 0.167277 0.197964 0.236462 0.292341 0.24074 0.303223 0.186687 0.216004 0.180331 0.193295 0.220631 0.163627 0.231513 0.227872 0.174596 0.438431 0.232857 0.173434 0.24075 0.205884 0.229944 0.21967 0.245279 0.216871 0.20197 0.192713 0.177422 0.296578 0.196594 0.185809 0.192957 0.201057 0.169487 0.27553 0.154566 0.288822 0.19408 0.295845 0.302895 0.168368 0.179606 0.307357 0.17545 0.209677 0.236071 0.232955 0.325679 0.169974 0.216848 0.189689 0.208012 0.199646 0.183753 0.220847 0.220814 0.147351 0.264607 0.202196 0.2079 0.298693 0.163476 0.143577 0.232592 0.186625 0.217039 0.185612 0.221597 0.347086 0.227663 0.173844 0.194194 0.18087 0.24537 0.338116 0.209235 0.192703 0.318752 0.293855 0.266708 0.214413 0.168105 0.249316 0.196055 0.192453 0.294101 0.265581 0.199001 0.273785 0.24714 0.267011 0.227785 0.192625 0.198373 0.218008 0.235256 0.280162 0.133567 0.140613 0.237799 0.257083 0.189969 0.251524 0.210618 0.205398 0.17587 0.180204 0.137759 0.183791 0.17364 0.167012 0.159311 0.256271 0.197267 0.128986 0.172898 0.146892 0.16842 0.213636 0.233459 0.189815 0.200319 0.192298 0.281148 0.294907 0.196121 0.213615 0.305068 0.205 0.185543 0.209811 0.222343 0.153974 0.220982 0.219983 0.179293 0.249573 0.351801 0.230316 0.275205 0.150744 0.231454 0.200757 0.301415 0.146716 0.357378 0.169709 0.219325 0.167663 0.365673 0.136454 0.272886 0.194938 0.163142 0.216449 0.145376 0.202919 0.237984 0.248437 0.145195 0.153015 0.193225 0.212384 0.21673 0.195015 0.188065 0.184078 0.226422 0.263999 0.194413 0.182799 0.18706 0.347263 0.239701 0.24118 0.268742 0.236538 0.224908 0.201437 0.242414 0.206753 0.2619 0.139622 0.209244 0.23448 0.191481 0.127425 0.261202 0.346355 0.229589 0.207934 0.248634 0.159952 0.24638 0.179538 0.219466 0.16505 0.207276 0.174766 0.237008 0.246077 0.167704 0.164729 0.203961 0.174964 0.194958 0.215422 0.160266 0.188479 0.197444 0.272562 0.203575 0.2259 0.277204 0.206028 0.178479 0.211094 0.319673 0.2113 0.223769 0.293053 0.208586 0.199428 0.225612 0.264203 0.152663 0.21918 0.290789 0.235908 0.262784 0.235217 0.157324 0.151208 0.182928 0.296855 0.186147 0.159248 0.239297 0.241843 0.340856 0.146117 0.203663 0.329354 0.249557 0.197849 0.229585 0.223439 0.20413 0.229131 0.171951 0.228079 0.157295 0.199271 0.271587 0.203292 0.281969 0.2444 0.230039 0.321699 0.165935 0.364062 0.252401 0.297667 0.184013 0.236464 0.204243 0.164873 0.259218 0.393396 0.211955 0.248085 0.160084 0.236716 0.15983 0.252671 0.239405 0.188662 0.264606 0.319808 0.307493 0.196324 0.181518 0.254667 0.264809 0.223117 0.26777 0.180863 0.169455 0.224828 0.207817 0.140732 0.213399 0.324769 0.191137 0.164256 0.217114 0.214011 0.32828 0.238256 0.143552 0.169232 0.186003 0.229799 0.261061 0.19808 0.239963 0.272574 0.207236 0.161954 0.177548 0.193625 0.165596 0.179549 0.214562 0.167227 0.151187 0.209401 0.241913 0.273214 0.218967 0.235485 0.175937 0.207597 0.185031 0.219831 0.187479 0.221324 0.306614 0.231653 0.19831 0.24108 0.328976 0.282177 0.158874 0.216567 0.229259 0.277304 0.179464 0.207069 0.14113 0.226342 0.145619 0.293699 0.248246 0.175565 0.272901 0.207309 0.387656 0.435315 0.216839 0.245956 0.239582 0.154019 0.260458 0.195431 0.149823 0.293539 0.210647 0.214458 0.358645 0.118122 0.275874 0.227726 0.176452 0.256971 0.151 0.214041 0.304766 0.224636 0.249624 0.166835 0.181667 0.307696 0.159636 0.169763 0.204124 0.219669 0.211665 0.207583 0.30489 0.217712 0.276536 0.150663 0.211808 0.203405 0.192382 0.235546 0.188585 0.176239 0.201623 0.304667 0.254138 0.191456 0.225764 0.233401 0.376006 0.344368 0.288585 0.180343 0.333063 0.178402 0.162344 0.227907 0.268659 0.232718 0.204081 0.206296 0.209226 0.192916 0.25015 0.205142 0.215824 0.208245 0.263261 0.171421 0.301258 0.212948 0.219135 0.159329 0.220561 0.216855 0.169559 0.172179 0.228201 0.233431 0.179697 0.276985 0.291083 0.165513 0.191762 0.206852 0.251609 0.150542 0.242458 0.196936 0.224343 0.21223 0.145229 0.186454 0.200116 0.212882 0.241305 0.112505 0.195894 0.236992 0.192106 0.192524 0.177171 0.348234 0.204171 0.204994 0.265007 0.291601 0.193835 0.177148 0.240226 0.209017 0.311512 0.215141 0.179632 0.248602 0.226228 0.201696 0.220132 0.252229 0.175802 0.203117 0.23226 0.168819 0.203104 0.155696 0.240561 0.274066 0.215836 0.15489 0.236 0.234238 0.217686 0.23863 0.206179 0.249432 0.193886 0.253006 0.198234 0.219076 0.254798 0.181505 0.206856 0.32696 0.166742 0.200603 0.200971 0.246308 0.303292 0.249063 0.173019 0.24702 0.227231 0.169647 0.212904 0.254934 0.255037 0.202556 0.148191 0.221891 0.222612 0.140609 0.202277 0.197141 0.21284 0.167761 0.163184 0.210675 0.227038 0.229011 0.298196 0.16029 0.211503 0.25007 0.277543 0.189126 0.239957 0.316071 0.187054 0.193545 0.190876 0.223157 0.365322 0.234077 0.296911 0.201239 0.235502 0.244947 0.225992 0.217926 0.322573 0.179331 0.263374 0.199866 0.180312 0.195591 0.294465 0.237645 0.202145 0.242913 0.210386 0.263999 0.204056 0.18286 0.203443 0.2247 0.167352 0.196426 0.202831 0.188345 0.18468 0.214904 0.22798 0.468077 0.153338 0.191507 0.185125 0.185146 0.201889 0.317922 0.192915 0.207596 0.283423 0.190196 0.285537 0.160485 0.129616 0.212227 0.202214 0.195337 0.176046 0.24785 0.146365 0.190598 0.205465 0.22706 0.150865 0.279961 0.204864 0.219338 0.245924 0.270775 0.143985 0.202874 0.237628 0.291821 0.245804 0.216056 0.285604 0.239524 0.218997 0.186178 0.259695 0.195732 0.224823 0.229557 0.18173 0.160023 0.260176 0.245332 0.194832 0.197061 0.193949 0.157516 0.165522 0.219592 0.331629 0.199656 0.322394 0.234908 0.21827 0.168228 0.267048 0.21352 0.270904 0.305966 0.239358 0.177989 0.18114 0.219485 0.212878 0.227175 0.299486 0.224959 0.241449 0.24692 0.385112 0.172421 0.166314 0.279355 0.202378 0.136757 0.191892 0.294383 0.174865 0.266785 0.223335 0.190861 0.24352 0.234891 0.244276 0.175541 0.222146 0.19296 0.21566 0.246103 0.293633 0.229102 0.33725 0.253867 0.218075 0.194736 0.183108 0.190547 0.140548 0.186624 0.239606 0.358013 0.25935 0.232492 0.245807 0.193667 0.161122 0.176418 0.237714 0.237863 0.232436 0.189605 0.221384 0.214452 0.178631 0.206126 0.249589 0.161993 0.170781 0.126471 0.261849 0.186842 0.210594 0.175542 0.225692 0.197006 0.286674 0.236582 0.195163 0.204741 0.271136 0.18414 0.232505 0.148779 0.146144 0.219079 0.228868 0.196311 0.244159 0.151122 0.18421 0.293976 0.19388 0.270681 0.191992 0.19549 0.174916 0.204016 0.232069 0.213043 0.240065 0.300391 0.176046 0.161706 0.189649 0.179124 0.202874 0.467081 0.22059 0.250633 0.140132 0.316238 0.253232 0.175695 0.218285 0.220841 0.209904 0.255766 0.163164 0.168085 0.178027 0.310219 0.269762 0.151308 0.309283 0.173199 0.171594 0.245815 0.164323 0.182936 0.191682 0.189303 0.15426 0.216314 0.379914 0.165325 0.16794 0.151221 0.244495 0.214393 0.246596 0.149193 0.275516 0.268829 0.137548 0.248692 0.180617 0.223001 0.265654 0.27262 0.16907 0.311469 0.418045 0.15488 0.194875 0.194258 0.169991 0.326691 0.155683 0.184986 0.208716 0.192315 0.299243 0.200334 0.149394 0.203992 0.239532 0.205563 0.321858 0.193249 0.20132 0.239567 0.156825 0.2774 0.381787 0.35785 0.202082 0.226624 0.286637 0.244074 0.161234 0.390729 0.216729 0.212425 0.199894 0.2207 0.281001 0.24932 0.273728 0.181905 0.215988 0.205149 0.207682 0.218455 0.241385 0.205827 0.184687 0.213545 0.217172 0.211174 0.228562 0.233369 0.157731 0.207448 0.16423 0.330905 0.180859 0.172804 0.397013 0.184331 0.279747 0.246054 0.201484 0.221644 0.3244 0.151236 0.161409 0.174604 0.203818 0.146524 0.170056 0.284361 0.24666 0.280192 0.329208 0.181953 0.163662 0.210672 0.205262 0.213747 0.182706 0.212833 0.229581 0.184294 0.306757 0.174432 0.293406 0.325316 0.191742 0.192447 0.221696 0.247693 0.165133 0.202366 0.15833 0.306401 0.209143 0.280897 0.238697 0.213136 0.187751 0.207594 0.211893 0.176591 0.22477 0.220267 0.197026 0.275449 0.197157 0.182913 0.217671 0.194453 0.241188 0.270775 0.262589 0.192895 0.176217 0.199603 0.182143 0.206158 0.207743 0.260882 0.255894 0.271596 0.238058 0.194117 0.128272 0.19302 0.207122 0.174719 0.252215 0.211388 0.161171 0.244758 0.179092 0.237612 0.187155 0.269728 0.149482 0.281395 0.185663 0.246065 0.287716 0.242477 0.166184 0.202454 0.213554 0.208457 0.210552 0.14398 0.215161 0.132728 0.188748 0.312969 0.347812 0.203051 0.338752 0.22179 0.168881 0.297266 0.329835 0.167859 0.287926 0.157087 0.251117 0.156516 0.142186 0.149646 0.216267 0.32152 0.329525 0.175946 0.182614 0.262386 0.216727 0.150641 0.285741 0.233647 0.225183 0.419613 0.179665 0.229122 0.225176 0.184388 0.246496 0.188243 0.17085 0.27503 0.19297 0.346195 0.188833 0.224257 0.212594 0.252832 0.177388 0.175431 0.271915 0.189173 0.229218 0.177333 0.138975 0.19347 0.170294 0.251635 0.208572 0.214356 0.2957 0.189839 0.247092 0.252045 0.288753 0.191632 0.227812 0.21789 0.281462 0.321454 0.35727 0.284309 0.201295 0.187946 0.21783 0.251502 0.18459 0.223785 0.231769 0.159275 0.255174 0.18606 0.174618 0.204641 0.170573 0.249083 0.187559 0.305846 0.191475 0.275309 0.242654 0.165841 0.169251 0.205466 0.229838 0.170291 0.222572 0.228612 0.236038 0.199379 0.201156 0.124251 0.275232 0.198599 0.194342 0.197445 0.189566 0.18129 0.218961 0.19765 0.243592 0.219266 0.197101 0.201392 0.19582 0.183541 0.216298 0.280882 0.18382 0.248951 0.208552 0.192414 0.158686 0.182694 0.192454 0.500638 0.213405 0.208632 0.199619 0.190697 0.214248 0.285199 0.140885 0.314417 0.280133 0.462125 0.231954 0.256618 0.31329 0.181585 0.31146 0.256131 0.186018 0.320828 0.18501 0.156544 0.20525 0.229831 0.196044 0.223575 0.168735 0.170017 0.241054 0.189162 0.23085 0.2255 0.407562 0.20615 0.195688 0.198861 0.173269 0.239074 0.429252 0.266407 0.152111 0.197567 0.43882 0.187293 0.198206 0.228058 0.148141 0.279034 0.173604 0.260208 0.192008 0.229907 0.197654 0.199112 0.1713 0.227572 0.215131 0.19613 0.214747 0.217658 0.238185 0.167865 0.22206 0.201318 0.194689 0.217341 0.186339 0.170162 0.213222 0.246195 0.266212 0.190687 0.196397 0.272658 0.19449 0.281267 0.181318 0.212607 0.25536 0.191718 0.316417 0.214117 0.172068 0.185053 0.149907 0.185545 0.191673 0.305081 0.161486 0.377373 0.238719 0.272344 0.239271 0.381949 0.183543 0.18502 0.123427 0.19392 0.176942 0.22449 0.20167 0.233519 0.135379 0.15717 0.174164 0.209512 0.252807 0.165604 0.215775 0.249171 0.252691 0.198674 0.209682 0.198446 0.194839 0.202653 0.184804 0.163108 0.213226 0.178571 0.19483 0.211225 0.20158 0.237051 0.168283 0.122038 0.245799 0.17989 0.266376 0.231679 0.233728 0.180929 0.176572 0.158854 0.256633 0.27149 0.163834 0.271086 0.204349 0.226359 0.193364 0.334167 0.201765 0.150435 0.259195 0.257097 0.309935 0.231378 0.341325 0.24328 0.251285 0.171986 0.179111 0.208423 0.212607 0.200593 0.384393 0.248546 0.294256 0.199687 0.150528 0.119705 0.139998 0.267272 0.283252 0.214955 0.283416 0.184194 0.19909 0.225762 0.225249 0.321455 0.329126 0.249132 0.328894 0.242243 0.193488 0.219321 0.166756 0.14049 0.142253 0.25548 0.199934 0.22337 0.199383 0.214055 0.250612 0.195737 0.176158 0.340546 0.150081 0.238243 0.202221 0.230653 0.199163 0.238826 0.195181 0.21218 0.188119 0.213793 0.221917 0.198496 0.171523 0.208021 0.203689 0.27087 0.197753 0.211966 0.166557 0.126906 0.192389 0.228077 0.165058 0.158951 0.206851 0.190268 0.258309 0.244561 0.252675 0.195132 0.185078 0.208613 0.244756 0.150271 0.281455 0.178624 0.11475 0.191 0.137471 0.205261 0.307997 0.273081 0.176642 0.154265 0.12617 0.215881 0.199388 0.225935 0.229192 0.279607 0.175209 0.309655 0.167447 0.266314 0.273724 0.235599 0.253238 0.368373 0.157468 0.202264 0.188677 0.216937 0.190687 0.254388 0.154962 0.156201 0.261991 0.214977 0.211906 0.221499 0.176977 0.196065 0.210634 0.219267 0.131424 0.151953 0.251483 0.27888 0.265352 0.192209 0.175631 0.274721 0.232022 0.177012 0.278178 0.216713 0.188588 0.184963 0.2058 0.269467 0.210208 0.230099 0.183752 0.226878 0.204455 0.183174 0.149866 0.245358 0.211169 0.205474 0.248674 0.176269 0.352567 0.18347 0.179079 0.161548 0.236859 0.375649 0.193158 0.226738 0.161634 0.210174 0.207757 0.221759 0.208394 0.223655 0.153863 0.31704 0.159847 0.209248 0.137941 0.133745 0.233025 0.33011 0.273993 0.226517 0.218394 0.291123 0.17579 0.256482 0.19338 0.248155 0.330072 0.255948 0.161672 0.217567 0.173189 0.212424 0.182322 0.326445 0.146103 0.267153 0.20662 0.194824 0.252306 0.228541 0.280331 0.206944 0.234907 0.16963 0.172965 0.190821 0.25979 0.246832 0.182695 0.222039 0.238074 0.25007 0.147418 0.18005 0.165103 0.190698 0.198994 0.199398 0.247777 0.17382 0.191645 0.192647 0.288628 0.175413 0.213123 0.230816 0.24668 0.248011 0.198778 0.291207 0.218072 0.185762 0.177643 0.263106 0.190104 0.176879 0.165668 0.327218 0.174 0.210574 0.289305 0.270406 0.20134 0.176822 0.26584 0.239028 0.201153 0.193388 0.194126 0.165557 0.217183 0.226376 0.167405 0.180808 0.131812 0.229994 0.212507 0.219203 0.169808 0.258754 0.237474 0.18592 0.290339 0.189381 0.183995 0.229041 0.217465 0.167814 0.248442 0.247976 0.221852 0.128478 0.243264 0.253365 0.259412 0.239568 0.210142 0.216653 0.240842 0.355557 0.247096 0.175604 0.217967 0.219845 0.181635 0.251259 0.183686 0.255949 0.207177 0.196353 0.239494 0.270308 0.19293 0.349375 0.209336 0.212577 0.179049 0.240141 0.190106 0.264814 0.235299 0.213348 0.230595 0.169097 0.21948 0.174257 0.194984 0.232705 0.168621 0.264862 0.309983 0.218021 0.217633 0.147346 0.121855 0.319091 0.17339 0.174501 0.190412 0.173877 0.204399 0.147158 0.264475 0.155916 0.211956 0.167467 0.269683 0.22363 0.249344 0.306502 0.184522 0.274235 0.240533 0.185065 0.218717 0.205016 0.279066 0.183254 0.220377 0.267355 0.249094 0.227858 0.237982 0.221577 0.232634 0.201605 0.19726 0.215769 0.216953 0.212984 0.436419 0.164546 0.196933 0.21806 0.184191 0.262096 0.21357 0.138956 0.176953 0.227409 0.285682 0.232094 0.174932 0.187615 0.398507 0.217882 0.231292 0.201083 0.205677 0.183435 0.301578 0.201618 0.197304 0.243681 0.188647 0.23247 0.214716 0.173421 0.18985 0.207925 0.247113 0.190513 0.186653 0.246466 0.202567 0.20121 0.199777 0.246056 0.373502 0.279706 0.205145 0.233426 0.280498 0.165361 0.20283 0.174212 0.197799 0.178052 0.207262 0.201194 0.254424 0.219293 0.256676 0.196115 0.184153 0.382053 0.26191 0.252188 0.285493 0.231394 0.252544 0.203045 0.203336 0.256068 0.168543 0.181576 0.330739 0.231885 0.159401 0.220153 0.31073 0.118172 0.21819 0.209155 0.239541 0.270095 0.198766 0.230958 0.286638 0.212382 0.250406 0.185916 0.254383 0.246141 0.181796 0.383203 0.208777 0.274149 0.195001 0.37781 0.272524 0.223798 0.140125 0.157907 0.274631 0.205127 0.302331 0.171109 0.294146 0.225321 0.20212 0.219387 0.21788 0.152611 0.177671 0.184441 0.179853 0.298839 0.224097 0.230564 0.194005 0.259329 0.231869 0.224748 0.198806 0.133686 0.172344 0.216872 0.185862 0.199504 0.158036 0.193558 0.198703 0.170047 0.161203 0.230348 0.21867 0.212289 0.187696 0.21699 0.160516 0.18228 0.266341 0.125623 0.179951 0.205344 0.204165 0.228807 0.128188 0.245349 0.254967 0.270173 0.28886 0.193988 0.155093 0.290844 0.159771 0.285379 0.173795 0.245165 0.147773 0.192166 0.260765 0.24263 0.19667 0.202291 0.142819 0.262632 0.206556 0.148579 0.199248 0.248227 0.255321 0.162296 0.160754 0.161587 0.189806 0.36664 0.279332 0.168264 0.192242 0.256311 0.247275 0.184111 0.18145 0.278315 0.123525 0.264556 0.238696 0.211341 0.18368 0.198309 0.221969 0.287223 0.132213 0.213208 0.235717 0.274711 0.161803 0.150197 0.20349 0.278816 0.216664 0.195247 0.244589 0.221768 0.181804 0.41886 0.242947 0.357496 0.283934 0.212022 0.285695 0.166787 0.217548 0.168436 0.193131 0.217258 0.275869 0.164649 0.195413 0.243083 0.277818 0.209279 0.231663 0.209859 0.176593 0.192832 0.229599 0.202668 0.237589 0.232839 0.165039 0.242642 0.181013 0.341582 0.25936 0.267922 0.230697 0.199647 0.413432 0.25759 0.422126 0.201545 0.327335 0.230144 0.18736 0.293437 0.276062 0.238636 0.150219 0.221126 0.196764 0.271353 0.208878 0.225424 0.241237 0.225814 0.167508 0.241214 0.18377 0.165063 0.183144 0.243671 0.398371 0.231227 0.222275 0.174135 0.209084 0.233265 0.158283 0.166834 0.214235 0.202569 0.275845 0.15292 0.204688 0.185051 0.278346 0.198631 0.208549 0.210011 0.196853 0.178872 0.19546 0.156793 0.209369 0.105762 0.189194 0.201094 0.16214 0.17899 0.185771 0.224502 0.263765 0.295179 0.163918 0.144799 0.150707 0.22468 0.26864 0.159725 0.206284 0.184721 0.173703 0.261767 0.232483 0.200421 0.151158 0.280614 0.214513 0.233102 0.262558 0.192175 0.229255 0.308165 0.191635 0.256869 0.204287 0.229947 0.33144 0.246567 0.138734 0.303031 0.300965 0.159272 0.186731 0.212411 0.227058 0.300833 0.281181 0.181054 0.262359 0.119216 0.247512 0.264299 0.184256 0.181758 0.303651 0.241389 0.320766 0.187947 0.233346 0.277409 0.234143 0.176839 0.225463 0.160394 0.194524 0.21968 0.24691 0.23866 0.171969 0.243747 0.224569 0.348882 0.219796 0.195765 0.300946 0.199066 0.193902 0.192068 0.388817 0.223463 0.215706 0.206895 0.187659 0.203483 0.1686 0.186164 0.289644 0.329109 0.202144 0.138972 0.301332 0.182837 0.256061 0.296326 0.249265 0.209701 0.170372 0.267045 0.162789 0.386598 0.271622 0.289973 0.268005 0.230502 0.156579 0.226197 0.204519 0.200671 0.187303 0.23213 0.23234 0.172755 0.186747 0.185705 0.155875 0.165421 0.166763 0.216772 0.223203 0.203208 0.184475 0.17437 0.192153 0.241058 0.186925 0.179969 0.159721 0.227982 0.252962 0.190131 0.195841 0.242469 0.17743 0.186279 0.201463 0.258776 0.193601 0.224904 0.260524 0.211174 0.345811 0.213611 0.3027 0.191521 0.163489 0.212334 0.202823 0.248435 0.199286 0.169113 0.163144 0.219038 0.270526 0.198079 0.160371 0.188559 0.210562 0.126418 0.185973 0.207851 0.235328 0.200776 0.234567 0.354781 0.14086 0.226509 0.182633 0.1707 0.26813 0.196998 0.19367 0.322938 0.180585 0.334027 0.176094 0.223819 0.227016 0.15708 0.355838 0.199307 0.188419 0.190273 0.220806 0.340899 0.236905 0.198136 0.19694 0.231235 0.276657 0.166744 0.18659 0.19478 0.260601 0.205281 0.174318 0.256651 0.218761 0.147146 0.143069 0.296025 0.254576 0.179746 0.180837 0.237495 0.158117 0.162576 0.20311 0.122281 0.220023 0.18955 0.139998 0.152963 0.289058 0.222232 0.222766 0.205125 0.209886 0.168051 0.244485 0.174826 0.163054 0.306316 0.219658 0.147174 0.186046 0.281012 0.17459 0.19709 0.243784 0.307569 0.214151 0.178029 0.314586 0.323819 0.174835 0.216216 0.150486 0.25321 0.186064 0.160567 0.243954 0.207025 0.169851 0.149235 0.173794 0.179439 0.23951 0.177618 0.350648 0.148429 0.178152 0.166672 0.225161 0.237496 0.222873 0.175038 0.184735 0.258875 0.130383 0.174851 0.248812 0.258573 0.234422 0.237588 0.212354 0.230822 0.164206 0.184652 0.176877 0.161302 0.130146 0.183324 0.283361 0.212937 0.241594 0.234739 0.208777 0.293253 0.198061 0.298024 0.195919 0.192039 0.2513 0.145163 0.200457 0.169142 0.303392 0.218553 0.168579 0.202667 0.213084 0.17036 0.184362 0.141618 0.223106 0.204943 0.170817 0.22133 0.171686 0.228034 0.359674 0.167 0.248886 0.182803 0.245203 0.236872 0.181689 0.212543 0.249907 0.156292 0.203628 0.151656 0.314555 0.177685 0.288319 0.312749 0.184744 0.223022 0.173226 0.223755 0.165383 0.174099 0.336876 0.215016 0.208422 0.345485 0.177131 0.189704 0.226849 0.176766 0.290566 0.204133 0.198504 0.222091 0.225141 0.189989 0.236368 0.175204 0.238492 0.161567 0.182682 0.305334 0.226717 0.213568 0.202532 0.172274 0.291138 0.200479 0.216679 0.194982 0.157246 0.214141 0.219505 0.197604 0.252118 0.228951 0.264628 0.209722 0.223879 0.190859 0.429809 0.143752 0.314365 0.283225 0.202865 0.183079 0.14333 0.181054 0.160161 0.215364 0.260416 0.193609 0.180964 0.180026 0.16354 0.220389 0.191486 0.166421 0.188078 0.13408 0.225343 0.257255 0.169971 0.271793 0.335194 0.288911 0.250154 0.228051 0.178087 0.157842 0.214961 0.138869 0.171814 0.217451 0.200927 0.258265 0.2045 0.184159 0.2009 0.258713 0.234502 0.286982 0.156707 0.144071 0.201821 0.22499 0.214637 0.211451 0.208804 0.233342 0.201494 0.21137 0.203667 0.198692 0.186745 0.323155 0.218859 0.264098 0.237679 0.242838 0.144071 0.216002 0.154463 0.205032 0.211065 0.287607 0.229187 0.200869 0.175263 0.21618 0.285803 0.18078 0.446582 0.215328 0.183307 0.275158 0.211631 0.197394 0.137332 0.160982 0.349922 0.163793 0.225294 0.237883 0.184748 0.172821 0.253585 0.294024 0.202352 0.204222 0.14024 0.169365 0.232889 0.175953 0.141385 0.258953 0.215406 0.143292 0.231925 0.196872 0.240916 0.11849 0.239344 0.251603 0.218708 0.291657 0.217763 0.154774 0.21608 0.19192 0.204587 0.24101 0.173281 0.249971 0.332573 0.200988 0.213825 0.286081 0.206289 0.194114 0.277931 0.248378 0.283997 0.199215 0.203312 0.226874 0.216773 0.23296 0.152109 0.235944 0.168171 0.248629 0.176229 0.254908 0.202514 0.220453 0.268119 0.332284 0.228824 0.231024 0.24188 0.179749 0.189311 0.186186 0.220214 0.242401 0.205381 0.219737 0.23302 0.179081 0.258617 0.174122 0.121867 0.294872 0.242165 0.195602 0.143097 0.197687 0.23496 0.170631 0.180277 0.165624 0.202823 0.18949 0.188559 0.204357 0.177874 0.170541 0.130961 0.250273 0.223708 0.134976 0.174945 0.195468 0.261209 0.117102 0.197626 0.197403 0.210465 0.34465 0.211623 0.189828 0.192123 0.209001 0.191215 0.21574 0.220139 0.318621 0.208756 0.170282 0.190123 0.252134 0.175163 0.186684 0.221696 0.190379 0.179558 0.216915 0.259437 0.206591 0.171042 0.501645 0.26996 0.195972 0.159227 0.163638 0.265251 0.157642 0.188959 0.190356 0.284554 0.306262 0.205798 0.254458 0.32179 0.242061 0.25434 0.246179 0.318914 0.203628 0.193663 0.225078 0.201086 0.182444 0.205206 0.170963 0.155689 0.21923 0.284781 0.209331 0.261622 0.315304 0.170665 0.124827 0.250503 0.220868 0.192542 0.192664 0.146765 0.146853 0.192706 0.134988 0.220909 0.254067 0.200703 0.197264 0.115836 0.337319 0.251373 0.206408 0.14324 0.312487 0.140528 0.283959 0.190666 0.215531 0.189622 0.169298 0.203786 0.201941 0.289374 0.135222 0.291114 0.190298 0.33031 0.190889 0.225563 0.361999 0.197424 0.180709 0.199083 0.200802 0.249315 0.19695 0.16729 0.349389 0.269186 0.191171 0.188841 0.257579 0.182657 0.157164 0.235888 0.171307 0.312696 0.222381 0.230178 0.226953 0.264559 0.190229 0.226454 0.290463 0.213227 0.276493 0.2058 0.217185 0.238557 0.198467 0.151959 0.264221 0.185462 0.234215 0.179905 0.162733 0.213678 0.271158 0.26317 0.170455 0.152885 0.189517 0.166008 0.240636 0.123727 0.264451 0.161999 0.143805 0.188298 0.12525 0.209591 0.183226 0.186141 0.181152 0.173764 0.158157 0.203808 0.238482 0.19667 0.142762 0.146719 0.241112 0.176952 0.25542 0.16379 0.154163 0.207015 0.200013 0.267304 0.231199 0.185044 0.151155 0.188749 0.230278 0.24422 0.236132 0.182759 0.168805 0.178461 0.332304 0.1883 0.23329 0.17471 0.237512 0.297814 0.211536 0.143329 0.273133 0.203012 0.143162 0.308723 0.181272 0.225193 0.209954 0.253244 0.20104 0.248438 0.280105 0.302327 0.241501 0.209462 0.187198 0.193359 0.181105 0.163995 0.221134 0.276325 0.174027 0.29009 0.120581 0.299282 0.274045 0.211164 0.198097 0.209936 0.22449 0.19924 0.197102 0.21402 0.184757 0.218041 0.212554 0.183618 0.258754 0.179852 0.21317 0.246505 0.250976 0.269804 0.213266 0.178009 0.17747 0.274657 0.209993 0.172023 0.238182 0.171864 0.219343 0.252844 0.264635 0.171186 0.227507 0.176221 0.233519 0.218263 0.195186 0.34878 0.282098 0.248399 0.200806 0.188406 0.220315 0.189711 0.275 0.214825 0.18602 0.206926 0.25621 0.233071 0.223869 0.23993 0.197311 0.231902 0.319789 0.185643 0.371651 0.216651 0.188854 0.165634 0.182272 0.200559 0.147457 0.282048 0.183172 0.159261 0.268221 0.197125 0.185207 0.224366 0.20628 0.217471 0.181227 0.234085 0.208541 0.2488 0.354997 0.237737 0.159965 0.223015 0.185665 0.274723 0.174021 0.224271 0.246009 0.233639 0.182445 0.328124 0.170606 0.210308 0.156281 0.178815 0.1644 0.153397 0.171319 0.188029 0.228266 0.178563 0.274043 0.212992 0.263464 0.145279 0.235497 0.263487 0.187364 0.215082 0.252416 0.205685 0.184029 0.217474 0.193864 0.241484 0.214969 0.229979 0.181683 0.216065 0.199281 0.208071 0.181978 0.26505 0.223725 0.245493 0.186266 0.239339 0.260986 0.179996 0.210973 0.17227 0.262604 0.241993 0.260217 0.19933 0.175278 0.1369 0.246127 0.189851 0.128652 0.233747 0.278713 0.192513 0.255965 0.245425 0.293855 0.205305 0.285832 0.142922 0.217816 0.218469 0.186685 0.227274 0.22935 0.203751 0.163666 0.223187 0.20752 0.196918 0.213304 0.158288 0.194809 0.219691 0.32449 0.171508 0.235497 0.179852 0.202007 0.241609 0.148492 0.300089 0.145702 0.211654 0.225143 0.236397 0.209042 0.225923 0.168575 0.236442 0.255477 0.240401 0.175028 0.176867 0.161772 0.16837 0.410136 0.252228 0.259688 0.209975 0.184025 0.385023 0.143365 0.199421 0.277868 0.171036 0.221617 0.179209 0.233005 0.176737 0.156535 0.303053 0.206151 0.178469 0.25301 0.155956 0.195283 0.262143 0.212707 0.302662 0.203276 0.219153 0.228539 0.300929 0.231991 0.214466 0.216566 0.25117 0.192572 0.19508 0.324824 0.240844 0.202024 0.17389 0.352978 0.354091 0.218055 0.204421 0.18949 0.252011 0.291228 0.210082 0.199527 0.193164 0.229342 0.224425 0.259376 0.238766 0.165486 0.21717 0.207105 0.264778 0.239814 0.156996 0.226565 0.30295 0.16224 0.216641 0.189413 0.228793 0.355266 0.282413 0.207294 0.249207 0.170236 0.174537 0.202768 0.212347 0.173465 0.292789 0.182058 0.31543 0.215545 0.194848 0.182341 0.155017 0.227711 0.2362 0.158428 0.162046 0.228195 0.30215 0.226786 0.234269 0.210335 0.30363 0.179749 0.186289 0.1584 0.284884 0.194456 0.283112 0.321784 0.190793 0.253713 0.188857 0.13456 0.154482 0.163592 0.255584 0.248724 0.186963 0.183638 0.203032 0.154564 0.260048 0.185893 0.213925 0.162361 0.255515 0.231686 0.202535 0.180769 0.227986 0.184995 0.222428 0.164618 0.259408 0.136741 0.311542 0.203577 0.310381 0.215054 0.241341 0.210574 0.288802 0.17491 0.226794 0.271013 0.22233 0.322467 0.269587 0.166406 0.233358 0.195375 0.163575 0.176662 0.133722 0.220357 0.212186 0.345567 0.242821 0.314387 0.16879 0.215123 0.196731 0.252172 0.222229 0.208491 0.290355 0.195732 0.2587 0.268071 0.185509 0.199031 0.230476 0.284888 0.193732 0.175201 0.157308 0.221176 0.214507 0.210285 0.165064 0.167519 0.205695 0.238124 0.190388 0.28537 0.213917 0.230683 0.295975 0.196343 0.372712 0.192025 0.175641 0.182771 0.176409 0.339944 0.206905 0.286988 0.131519 0.276372 0.307072 0.238285 0.274994 0.221441 0.246016 0.241662 0.20717 0.464896 0.348591 0.338836 0.228572 0.352842 0.205443 0.205608 0.170181 0.178639 0.214471 0.256027 0.270909 0.241985 0.222481 0.181267 0.190746 0.262209 0.200703 0.218821 0.194091 0.26145 0.285554 0.217124 0.114746 0.25884 0.203678 0.168161 0.090654 0.217256 0.21427 0.177119 0.202866 0.138309 0.191847 0.110786 0.259761 0.431887 0.203961 0.198078 0.122866 0.193064 0.259581 0.232543 0.19303 0.21411 0.329105 0.237285 0.234821 0.226121 0.223848 0.19881 0.181953 0.242494 0.215591 0.292657 0.184668 0.15741 0.238652 0.197309 0.157443 0.256777 0.20952 0.225065 0.297305 0.163356 0.205623 0.209582 0.181089 0.226529 0.201913 0.206438 0.260514 0.228858 0.273915 0.300757 0.261008 0.237729 0.275446 0.219549 0.183336 0.225388 0.196516 0.268603 0.333514 0.204554 0.22106 0.172634 0.198646 0.337056 0.192902 0.339463 0.275752 0.219014 0.245195 0.18664 0.215164 0.220002 0.196751 0.17245 0.181415 0.187697 0.303512 0.329403 0.140201 0.247604 0.243831 0.179762 0.222051 0.177775 0.244025 0.154792 0.221317 0.168746 0.312058 0.162696 0.217136 0.190701 0.175132 0.281535 0.226523 0.189851 0.273767 0.182704 0.135291 0.21826 0.226302 0.158032 0.154278 0.164089 0.21373 0.233244 0.170893 0.166625 0.199736 0.128845 0.29347 0.293732 0.176347 0.292482 0.170015 0.251576 0.349455 0.235361 0.190742 0.208484 0.21047 0.172552 0.231583 0.152673 0.209385 0.231826 0.220254 0.193807 0.20704 0.180929 0.276155 0.173324 0.197765 0.179259 0.299599 0.173409 0.149051 0.195033 0.206286 0.20448 0.194981 0.2483 0.251799 0.199401 0.176945 0.309454 0.199135 0.231153 0.335588 0.27305 0.163503 0.191324 0.240535 0.284556 0.172865 0.209949 0.347167 0.183086 0.273046 0.360472 0.196176 0.239448 0.167181 0.203187 0.174342 0.318126 0.216588 0.20717 0.248912 0.200331 0.187786 0.149565 0.187142 0.1926 0.172665 0.293655 0.173407 0.19523 0.167329 0.143534 0.143218 0.245918 0.218953 0.218182 0.205157 0.375845 0.271093 0.256663 0.210196 0.244144 0.259777 0.191483 0.176498 0.27928 0.21343 0.163703 0.228986 0.192467 0.36256 0.163857 0.26308 0.149679 0.298121 0.165695 0.179313 0.12572 0.228344 0.210502 0.235685 0.122338 0.203112 0.199908 0.219254 0.216279 0.25331 0.16518 0.395939 0.202227 0.174687 0.223462 0.296872 0.226468 0.356442 0.13428 0.182245 0.326965 0.180138 0.329666 0.251507 0.184896 0.333592 0.295372 0.159469 0.170919 0.279634 0.189241 0.170971 0.170095 0.158769 0.198 0.195033 0.190289 0.186088 0.223603 0.199568 0.172234 0.308944 0.222201 0.346218 0.271309 0.171844 0.199653 0.240633 0.171012 0.197335 0.282596 0.263961 0.200192 0.232574 0.195065 0.270053 0.160901 0.240141 0.175691 0.196508 0.187003 0.332809 0.240299 0.168927 0.169085 0.207579 0.179292 0.227267 0.299079 0.231591 0.176554 0.260521 0.20063 0.255057 0.208833 0.21171 0.278779 0.224352 0.178935 0.268661 0.16866 0.169278 0.263294 0.338458 0.209924 0.257786 0.182322 0.181774 0.228475 0.208987 0.254731 0.189326 0.182059 0.154651 0.195495 0.232516 0.159186 0.214167 0.165113 0.317482 0.18848 0.253855 0.300697 0.206124 0.244244 0.222253 0.205431 0.207785 0.309662 0.211591 0.245712 0.258039 0.257966 0.187507 0.257138 0.335682 0.248959 0.222009 0.210486 0.205642 0.121298 0.234114 0.241639 0.404417 0.179223 0.234147 0.313622 0.165898 0.298745 0.276124 0.203421 0.189765 0.2193 0.216875 0.27453 0.229858 0.209953 0.171266 0.225963 0.189871 0.233451 0.2491 0.214006 0.1828 0.196804 0.247383 0.223138 0.16321 0.154887 0.234167 0.144434 0.131321 0.211378 0.174887 0.254635 0.218064 0.215036 0.152961 0.289504 0.227856 0.204409 0.127185 0.231419 0.269333 0.26766 0.199221 0.199924 0.224476 0.287791 0.168668 0.192804 0.291066 0.253214 0.243709 0.342811 0.207761 0.165106 0.242878 0.171123 0.215536 0.280293 0.140182 0.15807 0.269492 0.181769 0.245917 0.164423 0.297937 0.204483 0.197784 0.27994 0.285018 0.215528 0.179153 0.234578 0.24483 0.190211 0.224971 0.206825 0.168653 0.134115 0.148291 0.12478 0.229451 0.197181 0.336516 0.20806 0.210325 0.185104 0.182957 0.327843 0.224949 0.252987 0.314963 0.28048 0.255817 0.177379 0.237125 0.219295 0.268752 0.189654 0.266634 0.206496 0.334092 0.209089 0.217758 0.205783 0.216384 0.259202 0.259926 0.264688 0.17067 0.208718 0.248145 0.143282 0.20054 0.165773 0.180215 0.183102 0.152299 0.243778 0.209836 0.246815 0.238274 0.161671 0.265871 0.251926 0.153215 0.228783 0.186661 0.167842 0.16771 0.170147 0.158816 0.302555 0.192232 0.192336 0.25727 0.237058 0.307732 0.190332 0.185044 0.142403 0.202639 0.187598 0.195093 0.259306 0.243402 0.182679 0.22427 0.17152 0.299703 0.222131 0.153913 0.245252 0.219103 0.189701 0.148879 0.277988 0.217713 0.160868 0.227225 0.213095 0.208212 0.260062 0.270319 0.259159 0.141911 0.223294 0.199436 0.252678 0.168884 0.243581 0.148901 0.218911 0.182236 0.242581 0.402014 0.160547 0.319185 0.213691 0.238096 0.105544 0.203344 0.248525 0.197623 0.477808 0.181971 0.332076 0.248332 0.181784 0.257149 0.255498 0.156827 0.209571 0.242092 0.207898 0.207192 0.26641 0.266295 0.191988 0.200299 0.339886 0.19164 0.209346 0.158006 0.210291 0.192424 0.188972 0.262327 0.195434 0.180501 0.199605 0.187758 0.208477 0.12671 0.235047 0.253069 0.161898 0.215644 0.233521 0.25063 0.194122 0.170941 0.224701 0.241997 0.250432 0.169224 0.202087 0.194273 0.206823 0.165189 0.168069 0.18903 0.195275 0.186122 0.235602 0.166501 0.252385 0.25749 0.208713 0.243288 0.247792 0.185648 0.159335 0.12413 0.22926 0.206055 0.20589 0.203368 0.19039 0.241994 0.19213 0.238306 0.240508 0.182489 0.225145 0.262421 0.209905 0.216248 0.540176 0.205791 0.24367 0.206964 0.213706 0.320627 0.219964 0.186603 0.174838 0.194477 0.22677 0.201604 0.216212 0.207629 0.174486 0.158072 0.168369 0.134456 0.287331 0.319547 0.193683 0.264836 0.286839 0.214981 0.238267 0.206152 0.188373 0.22934 0.182861 0.230555 0.206488 0.175409 0.308975 0.269566 0.236139 0.256193 0.172353 0.257699 0.161506 0.210746 0.317839 0.208974 0.160141 0.263174 0.207778 0.259005 0.167855 0.212164 0.253977 0.163654 0.30109 0.313113 0.460769 0.173919 0.206526 0.242427 0.163373 0.202115 0.24515 0.466291 0.187501 0.134004 0.22497 0.269431 0.208964 0.183101 0.151113 0.219385 0.244013 0.189846 0.21765 0.303815 0.183224 0.190052 0.254315 0.201537 0.205292 0.168518 0.197681 0.199098 0.24767 0.232973 0.380487 0.273205 0.168564 0.274222 0.206145 0.216881 0.255524 0.288727 0.16784 0.264114 0.182503 0.183611 0.193229 0.191744 0.118163 0.177087 0.26139 0.232637 0.250415 0.207291 0.295834 0.272788 0.303743 0.337968 0.191001 0.264071 0.200039 0.209347 0.199019 0.179353 0.154273 0.255297 0.177726 0.273901 0.211305 0.183046 0.197316 0.2507 0.193107 0.181948 0.310127 0.196596 0.197224 0.395398 0.215044 0.265951 0.154501 0.246751 0.210341 0.137854 0.310951 0.172055 0.27591 0.270509 0.191731 0.218035 0.176575 0.205419 0.221115 0.389419 0.189492 0.229624 0.182742 0.446837 0.247213 0.211791 0.23772 0.17336 0.346717 0.185928 0.197208 0.241445 0.355936 0.160938 0.265204 0.198165 0.275132 0.193143 0.208758 0.150197 0.163743 0.179509 0.154106 0.233457 0.200054 0.239811 0.137288 0.156569 0.211559 0.269386 0.28683 0.169548 0.201226 0.183666 0.205629 0.214825 0.193031 0.168249 0.170662 0.374541 0.259995 0.178047 0.218638 0.192267 0.21265 0.174216 0.270212 0.204 0.2672 0.219276 0.295567 0.20592 0.167768 0.221128 0.19305 0.339972 0.152905 0.273863 0.246552 0.124779 0.249199 0.199724 0.349083 0.215893 0.190107 0.254459 0.185824 0.275209 0.190787 0.191593 0.194549 0.257641 0.204058 0.228406 0.349713 0.233201 0.223586 0.175855 0.244742 0.176684 0.179859 0.157642 0.408409 0.173511 0.328566 0.16844 0.187665 0.171892 0.182528 0.273838 0.24322 0.254941 0.10193 0.241637 0.267092 0.188769 0.230194 0.183254 0.216224 0.16769 0.165674 0.183661 0.223891 0.239894 0.238875 0.155985 0.196363 0.238184 0.169539 0.216279 0.223838 0.166597 0.154246 0.265795 0.379901 0.209198 0.175335 0.264579 0.278035 0.26553 0.211159 0.188305 0.211952 0.224225 0.19257 0.240821 0.181257 0.159756 0.208476 0.205022 0.252466 0.226763 0.176519 0.271866 0.260536 0.259053 0.252391 0.180048 0.226113 0.146726 0.187295 0.158863 0.161958 0.241177 0.201265 0.251769 0.255701 0.234205 0.233738 0.260614 0.17527 0.232941 0.223398 0.176233 0.225006 0.177648 0.258186 0.191544 0.174483 0.278282 0.217071 0.1761 0.175901 0.209069 0.249429 0.217537 0.180193 0.19395 0.227226 0.266155 0.230678 0.194707 0.22616 0.256084 0.194781 0.220424 0.209931 0.156535 0.247741 0.397576 0.188399 0.256389 0.210158 0.225878 0.170463 0.206964 0.200825 0.186978 0.18895 0.235649 0.25681 0.194926 0.21992 0.290316 0.232602 0.399267 0.209874 0.198561 0.1703 0.235291 0.234203 0.19555 0.23768 0.123609 0.457772 0.201969 0.20727 0.231893 0.166245 0.188108 0.187867 0.294855 0.15453 0.189593 0.235589 0.186962 0.281624 0.23859 0.225428 0.131547 0.198658 0.145243 0.259392 0.199326 0.185031 0.189581 0.22068 0.204433 0.163089 0.252957 0.272254 0.280335 0.26106 0.328377 0.258515 0.31174 0.286876 0.270716 0.204102 0.200272 0.22126 0.172005 0.23926 0.243154 0.187853 0.208678 0.215039 0.167673 0.178291 0.263558 0.234307 0.213772 0.21136 0.181643 0.123811 0.190611 0.231055 0.149367 0.166406 0.176879 0.181729 0.319519 0.141069 0.195825 0.209817 0.155076 0.206572 0.246775 0.263615 0.165684 0.336758 0.306974 0.527594 0.250815 0.325497 0.170828 0.21048 0.203849 0.192787 0.195606 0.202423 0.279123 0.135195 0.342053 0.188905 0.369816 0.261846 0.209134 0.205255 0.224475 0.368175 0.244249 0.151175 0.353779 0.206212 0.177002 0.128975 0.222092 0.307529 0.214696 0.201308 0.161695 0.213077 0.232342 0.178645 0.198269 0.175379 0.216008 0.230821 0.157593 0.192836 0.160966 0.11904 0.245928 0.32566 0.157782 0.147065 0.394299 0.165256 0.154894 0.256378 0.213929 0.203334 0.196726 0.223107 0.259349 0.343934 0.227365 0.210928 0.193339 0.39211 0.171857 0.142311 0.209441 0.277429 0.213662 0.201923 0.19881 0.285861 0.188638 0.179379 0.165367 0.205448 0.198357 0.26013 0.230029 0.199217 0.167233 0.197541 0.286459 0.222701 0.229297 0.196118 0.246654 0.20813 0.188666 0.327729 0.355238 0.216375 0.22581 0.193435 0.237392 0.195793 0.160329 0.329073 0.172684 0.236808 0.184816 0.219824 0.211848 0.263635 0.162626 0.187984 0.293432 0.311979 0.3178 0.257047 0.237085 0.232011 0.239337 0.17892 0.299658 0.249615 0.207903 0.224048 0.229898 0.178619 0.15374 0.162224 0.175235 0.289913 0.187357 0.272723 0.176618 0.254771 0.221882 0.254459 0.23308 0.246895 0.175373 0.160594 0.178045 0.249088 0.149333 0.218251 0.211041 0.149786 0.142686 0.223219 0.194424 0.185197 0.244957 0.185858 0.192203 0.425563 0.188819 0.201741 0.172216 0.227748 0.2062 0.18766 0.160195 0.169916 0.158268 0.194729 0.226156 0.23307 0.275178 0.186611 0.172402 0.171729 0.14865 0.220234 0.182381 0.314873 0.270404 0.199285 0.173861 0.295708 0.228307 0.173495 0.140642 0.18188 0.236632 0.198261 0.128797 0.193671 0.161231 0.191875 0.201848 0.169082 0.285378 0.164152 0.199739 0.216397 0.20314 0.135085 0.377275 0.235958 0.233594 0.248969 0.21558 0.233071 0.296976 0.317663 0.199054 0.200904 0.203591 0.272712 0.20562 0.229844 0.198919 0.173559 0.245451 0.216544 0.268319 0.183147 0.17217 0.172302 0.406091 0.197955 0.209985 0.241399 0.19958 0.204133 0.189903 0.397273 0.251105 0.207756 0.237894 0.225179 0.205683 0.207361 0.196289 0.193549 0.24935 0.182491 0.226166 0.169988 0.23275 0.248302 0.25724 0.273759 0.165377 0.098915 0.338216 0.223166 0.234454 0.18924 0.268573 0.212615 0.179511 0.19068 0.169924 0.182105 0.15406 0.221247 0.208914 0.1582 0.178257 0.275024 0.312395 0.134075 0.244115 0.279374 0.370545 0.274126 0.242443 0.225217 0.277731 0.162947 0.223526 0.224717 0.248834 0.209298 0.207732 0.181654 0.16766 0.161945 0.213293 0.223429 0.121455 0.163773 0.185163 0.271695 0.232846 0.19429 0.223788 0.317242 0.208767 0.251174 0.176816 0.286296 0.274216 0.274929 0.143382 0.175697 0.140976 0.191919 0.260682 0.240579 0.162166 0.116668 0.22352 0.196033 0.180455 0.228389 0.251691 0.185757 0.238683 0.197089 0.25043 0.288489 0.195404 0.217182 0.186253 0.237515 0.339667 0.232391 0.191692 0.325655 0.25096 0.362071 0.177756 0.222784 0.278521 0.21037 0.220363 0.196196 0.245304 0.18137 0.174182 0.183766 0.243505 0.229452 0.250478 0.215309 0.203817 0.23809 0.168082 0.204936 0.342289 0.141487 0.157019 0.232776 0.15365 0.490528 0.310402 0.182766 0.249792 0.168763 0.161716 0.193928 0.211811 0.187246 0.267753 0.211879 0.335345 0.165742 0.24347 0.160575 0.258528 0.149945 0.191036 0.162569 0.197269 0.190389 0.136704 0.143561 0.220872 0.187402 0.223854 0.303487 0.369841 0.110557 0.219505 0.251509 0.224862 0.2839 0.155867 0.205537 0.190934 0.278065 0.218655 0.200107 0.225979 0.202332 0.184677 0.320542 0.274273 0.195115 0.144885 0.193609 0.464871 0.173903 0.216641 0.189533 0.142206 0.202147 0.233966 0.108095 0.187338 0.211921 0.364899 0.203387 0.277642 0.248064 0.291652 0.246227 0.178665 0.235936 0.156296 0.178041 0.184175 0.229676 0.199023 0.247537 0.157102 0.210861 0.225169 0.14876 0.251722 0.113004 0.212633 0.177998 0.207121 0.28401 0.212474 0.304439 0.236093 0.193174 0.169887 0.314368 0.264054 0.210679 0.213791 0.231223 0.218625 0.227426 0.187137 0.236187 0.196885 0.242011 0.250152 0.196026 0.278625 0.370556 0.201705 0.243548 0.287354 0.195456 0.226031 0.259651 0.184634 0.235266 0.19205 0.205542 0.215834 0.230542 0.183295 0.321022 0.192092 0.237502 0.179572 0.267025 0.141861 0.236472 0.251965 0.23284 0.202041 0.190051 0.195108 0.186199 0.251308 0.236449 0.216064 0.200725 0.171022 0.157262 0.242488 0.169859 0.252358 0.178334 0.20255 0.286001 0.193871 0.186979 0.202988 0.249979 0.240501 0.158223 0.230334 0.196199 0.152446 0.186987 0.170738 0.223024 0.185727 0.144443 0.161956 0.245531 0.233978 0.181626 0.358526 0.250063 0.239376 0.225559 0.232259 0.235513 0.269775 0.185061 0.21909 0.168888 0.173358 0.174292 0.229101 0.194399 0.162023 0.176895 0.318466 0.233834 0.279324 0.212827 0.205566 0.304255 0.174706 0.223403 0.224137 0.203915 0.207566 0.221968 0.088753 0.249842 0.193459 0.294399 0.182847 0.223208 0.096061 0.171453 0.282141 0.346551 0.258221 0.286992 0.242341 0.191335 0.196095 0.277995 0.220015 0.172512 0.211422 0.207289 0.202438 0.321829 0.379393 0.274962 0.140109 0.151798 0.187075 0.260063 0.374699 0.190313 0.29192 0.219183 0.17191 0.2105 0.169294 0.198815 0.21609 0.160268 0.194534 0.127866 0.189277 0.224656 0.185645 0.279591 0.182601 0.231823 0.186355 0.1593 0.201839 0.170545 0.145651 0.179534 0.277366 0.179803 0.195531 0.175733 0.179755 0.155647 0.24337 0.154459 0.198519 0.288546 0.190605 0.252553 0.207173 0.236407 0.206781 0.249331 0.254919 0.212318 0.411522 0.178968 0.277566 0.232302 0.275009 0.344095 0.273508 0.234298 0.201894 0.269594 0.209249 0.177391 0.214538 0.301977 0.227368 0.230508 0.192078 0.225773 0.217939 0.373135 0.285759 0.361868 0.185108 0.360451 0.378654 0.161185 0.22916 0.20864 0.273286 0.20559 0.268168 0.271503 0.192511 0.426803 0.210799 0.199576 0.193145 0.386552 0.137295 0.238656 0.192565 0.172085 0.210062 0.173358 0.289547 0.28505 0.199695 0.184381 0.188144 0.193228 0.283893 0.193614 0.169621 0.161452 0.133622 0.153812 0.265455 0.243528 0.303775 0.219735 0.165844 0.279149 0.186457 0.197193 0.221084 0.169306 0.240252 0.204652 0.226752 0.163234 0.212561 0.157244 0.27439 0.191572 0.163093 0.148775 0.265265 0.203641 0.182917 0.260711 0.182975 0.222801 0.300554 0.189627 0.184342 0.253835 0.150016 0.158248 0.191396 0.219233 0.215359 0.218033 0.164029 0.274937 0.162432 0.137591 0.220756 0.228723 0.341629 0.205955 0.215326 0.226441 0.293822 0.265003 0.181401 0.150407 0.165545 0.25203 0.169635 0.278177 0.369532 0.257106 0.221651 0.21259 0.20668 0.347824 0.227434 0.177951 0.193787 0.246063 0.196656 0.359207 0.292573 0.224628 0.20487 0.148671 0.262223 0.236483 0.243261 0.207202 0.142174 0.349264 0.237992 0.21676 0.185736 0.283342 0.274717 0.236056 0.187038 0.24598 0.294308 0.259239 0.176251 0.1467 0.268717 0.234459 0.225224 0.263322 0.227612 0.180675 0.216399 0.232389 0.172411 0.17518 0.23668 0.203241 0.139013 0.319762 0.154138 0.150715 0.207962 0.194017 0.169975 0.166947 0.377396 0.258774 0.227265 0.24291 0.249282 0.143905 0.166956 0.141315 0.208815 0.182842 0.155643 0.227787 0.178278 0.209375 0.142761 0.200766 0.190233 0.136419 0.222003 0.187577 0.193712 0.197546 0.253246 0.242911 0.184101 0.312742 0.277462 0.148381 0.246356 0.190379 0.188296 0.199285 0.318135 0.283878 0.16635 0.214862 0.133179 0.352773 0.186894 0.186877 0.267633 0.229096 0.234245 0.215277 0.249408 0.233102 0.218889 0.186419 0.204871 0.170585 0.264432 0.222996 0.203304 0.392731 0.198041 0.167763 0.240796 0.27973 0.15796 0.171424 0.32398 0.156718 0.186868 0.245371 0.204099 0.265666 0.20287 0.178307 0.222872 0.174184 0.249078 0.302785 0.276513 0.232932 0.201074 0.251648 0.180197 0.163417 0.204884 0.255073 0.202389 0.190785 0.229334 0.239452 0.320649 0.170772 0.210247 0.183844 0.266649 0.153158 0.212805 0.217256 0.167461 0.253943 0.158617 0.273762 0.318889 0.167473 0.262071 0.242275 0.198036 0.234503 0.137775 0.13594 0.242276 0.282757 0.347099 0.188057 0.189146 0.140719 0.223594 0.217456 0.217284 0.225935 0.270769 0.179221 0.199647 0.238061 0.197643 0.198492 0.209936 0.25696 0.256045 0.180395 0.223269 0.116338 0.197279 0.184134 0.227056 0.237045 0.292633 0.219635 0.294891 0.290629 0.180404 0.195253 0.291072 0.219438 0.375064 0.155547 0.262727 0.247035 0.327856 0.185444 0.196003 0.281036 0.162528 0.13747 0.358555 0.200132 0.175509 0.22867 0.150176 0.412754 0.171253 0.249128 0.232318 0.175373 0.22581 0.175541 0.219732 0.19093 0.269115 0.141784 0.173006 0.171807 0.244987 0.1714 0.222512 0.195259 0.192813 0.189153 0.245387 0.230273 0.154946 0.200435 0.173859 0.169449 0.226664 0.211993 0.181231 0.294424 0.212314 0.219114 0.198305 0.20975 0.222294 0.196926 0.225775 0.176704 0.27994 0.224388 0.306983 0.198843 0.15442 0.171955 0.143082 0.162974 0.176326 0.188673 0.215167 0.15992 0.270948 0.463678 0.28495 0.272151 0.280968 0.181386 0.15657 0.202419 0.206873 0.175462 0.195563 0.157906 0.206897 0.267344 0.212955 0.194211 0.204871 0.280065 0.214081 0.244374 0.207232 0.202857 0.186614 0.219421 0.235739 0.431553 0.208895 0.202955 0.184321 0.186393 0.259844 0.264475 0.252561 0.306136 0.202955 0.248887 0.177303 0.21274 0.21066 0.1694 0.354717 0.232522 0.198079 0.248241 0.148794 0.238538 0.110873 0.227005 0.181392 0.207083 0.233629 0.252301 0.205346 0.180411 0.212446 0.171949 0.274921 0.251163 0.193537 0.242332 0.169294 0.233582 0.20386 0.17522 0.223019 0.188129 0.168021 0.216456 0.193687 0.19712 0.18819 0.171763 0.272296 0.183216 0.248951 0.209132 0.246872 0.122646 0.200693 0.249302 0.221751 0.242066 0.15814 0.219436 0.228066 0.155761 0.182741 0.361286 0.206905 0.291544 0.220543 0.298073 0.180134 0.20499 0.175083 0.240466 0.183995 0.178101 0.299695 0.120743 0.262437 0.117968 0.177341 0.300944 0.264628 0.293538 0.20137 0.244379 0.219139 0.178256 0.246544 0.218049 0.217753 0.260291 0.211194 0.330399 0.168995 0.264492 0.180977 0.243864 0.163465 0.276191 0.184281 0.127348 0.2158 0.359457 0.216347 0.121725 0.154752 0.262249 0.236761 0.222906 0.294766 0.248533 0.18295 0.25606 0.179794 0.169834 0.268469 0.245199 0.206827 0.20919 0.150917 0.197784 0.2499 0.222284 0.158652 0.219633 0.252849 0.227636 0.192779 0.32751 0.233469 0.188836 0.180226 0.176018 0.206599 0.218757 0.191453 0.214423 0.186228 0.159221 0.173291 0.211668 0.204753 0.193347 0.257879 0.267204 0.15694 0.180571 0.21127 0.344716 0.190603 0.186448 0.220242 0.192292 0.307377 0.230691 0.362529 0.165496 0.29838 0.220517 0.20204 0.179018 0.243662 0.190768 0.211271 0.135845 0.197908 0.169435 0.251787 0.243309 0.255031 0.177165 0.203641 0.229778 0.211269 0.186189 0.203743 0.163038 0.284578 0.268213 0.297368 0.171112 0.179603 0.20522 0.234351 0.215992 0.176511 0.243138 0.266898 0.236454 0.206651 0.187719 0.207248 0.356639 0.244843 0.151586 0.194334 0.147957 0.167799 0.265497 0.424203 0.239039 0.202459 0.249663 0.208696 0.281424 0.141354 0.18752 0.181494 0.20147 0.202709 0.262969 0.271244 0.187395 0.209409 0.185041 0.14385 0.258888 0.12279 0.295727 0.132667 0.226017 0.226318 0.195631 0.285143 0.185839 0.225775 0.248486 0.176894 0.173348 0.151717 0.212228 0.163514 0.181638 0.281313 0.345527 0.107474 0.200892 0.407568 0.35165 0.191394 0.293535 0.184583 0.165158 0.231882 0.210556 0.420498 0.260058 0.221328 0.281724 0.190141 0.202823 0.153759 0.223966 0.247141 0.197707 0.173407 0.203382 0.153306 0.209584 0.196085 0.167172 0.208077 0.254082 0.247712 0.182263 0.149306 0.230919 0.287038 0.302203 0.180506 0.216541 0.262114 0.320172 0.21168 0.141353 0.247781 0.299682 0.215815 0.233467 0.322161 0.229986 0.18542 0.228325 0.197678 0.207226 0.332556 0.265135 0.223999 0.188669 0.236875 0.183116 0.18977 0.169895 0.154508 0.189371 0.231907 0.187178 0.200621 0.253981 0.225104 0.236709 0.205909 0.315676 0.259248 0.197641 0.22737 0.208982 0.248245 0.178814 0.234381 0.207288 0.19174 0.232287 0.177414 0.235915 0.153784 0.198885 0.190946 0.229621 0.312623 0.263937 0.221852 0.342382 0.220543 0.288767 0.240511 0.186943 0.206139 0.280565 0.242417 0.227555 0.145336 0.165004 0.231888 0.175812 0.408254 0.194806 0.186429 0.301724 0.320194 0.226117 0.272572 0.233934 0.317806 0.378519 0.181388 0.228848 0.231133 0.263665 0.221241 0.136513 0.201708 0.291472 0.218989 0.27677 0.181385 0.169237 0.189938 0.214726 0.175475 0.16836 0.179145 0.293327 0.246303 0.147295 0.203067 0.193309 0.208381 0.281613 0.255552 0.176587 0.270092 0.27553 0.213131 0.158483 0.256466 0.155397 0.254824 0.143038 0.250233 0.222588 0.182444 0.253569 0.175405 0.202535 0.174837 0.216567 0.297754 0.260229 0.248624 0.247194 0.180623 0.172371 0.146556 0.174898 0.184779 0.203909 0.251603 0.21011 0.202728 0.170948 0.190342 0.194354 0.187405 0.306627 0.267718 0.209022 0.241381 0.284768 0.204317 0.184429 0.254437 0.236021 0.074317 0.181838 0.19945 0.233907 0.200717 0.363747 0.208368 0.166544 0.166955 0.239156 0.152733 0.336224 0.118743 0.222889 0.1529 0.206766 0.179998 0.191996 0.342853 0.22529 0.218094 0.185094 0.192092 0.183332 0.210786 0.36687 0.268834 0.160859 0.236025 0.198284 0.313206 0.196555 0.219606 0.170813 0.198293 0.273084 0.19466 0.159703 0.246996 0.190374 0.198494 0.206244 0.365195 0.181286 0.188112 0.305825 0.236418 0.225962 0.192981 0.263834 0.297111 0.189317 0.242031 0.334061 0.197572 0.1993 0.23721 0.211952 0.184887 0.229038 0.260412 0.28913 0.189259 0.207824 0.164617 0.195425 0.171096 0.152936 0.220523 0.24777 0.197411 0.196559 0.375603 0.167178 0.196064 0.292051 0.225195 0.251369 0.136566 0.231689 0.235507 0.249425 0.227507 0.285263 0.229392 0.387813 0.346661 0.201474 0.16521 0.190641 0.199971 0.222765 0.233811 0.238033 0.220849 0.235627 0.255644 0.125012 0.213936 0.208117 0.200869 0.191227 0.196656 0.20334 0.202116 0.17548 0.225106 0.168389 0.188027 0.151143 0.170093 0.160696 0.168961 0.314327 0.208723 0.19549 0.189336 0.197034 0.179828 0.190192 0.156065 0.29459 0.21329 0.201481 0.19949 0.165935 0.208978 0.212242 0.190685 0.1925 0.227672 0.296318 0.219368 0.240196 0.184237 0.21037 0.234581 0.162686 0.163835 0.216172 0.144086 0.157265 0.207748 0.244213 0.206886 0.259833 0.203372 0.146622 0.18474 0.206043 0.254989 0.205886 0.2311 0.20518 0.204631 0.179032 0.231315 0.165069 0.20522 0.25461 0.246356 0.315977 0.294263 0.247013 0.190559 0.26503 0.221692 0.219628 0.247216 0.189708 0.193886 0.131881 0.130507 0.22119 0.186963 0.194122 0.20259 0.209059 0.169252 0.176605 0.207396 0.341381 0.185341 0.322986 0.166199 0.203024 0.204377 0.130687 0.197629 0.162975 0.222793 0.185517 0.184988 0.357001 0.186238 0.109411 0.198444 0.159836 0.41296 0.122621 0.191908 0.153376 0.190775 0.196258 0.222234 0.163967 0.204817 0.169363 0.198789 0.175383 0.203516 0.181658 0.141918 0.287021 0.147793 0.260059 0.158976 0.214911 0.232359 0.253 0.168091 0.289509 0.19543 0.201236 0.231431 0.226792 0.346996 0.197675 0.274941 0.249914 0.164344 0.171302 0.211441 0.133062 0.293647 0.159759 0.300504 0.178863 0.238036 0.2525 0.269293 0.261202 0.174673 0.178768 0.302219 0.203547 0.163688 0.194707 0.25241 0.261717 0.200063 0.289652 0.218187 0.195931 0.248198 0.206047 0.178293 0.222482 0.171605 0.219572 0.297217 0.194187 0.181902 0.239709 0.165223 0.231119 0.187439 0.225086 0.22436 0.13249 0.202765 0.232888 0.233733 0.139695 0.237481 0.196912 0.162418 0.231433 0.216498 0.31263 0.181244 0.298514 0.179115 0.254686 0.26849 0.18709 0.1586 0.177437 0.173686 0.179008 0.25588 0.238098 0.241441 0.323741 0.149472 0.178366 0.165532 0.199258 0.198575 0.208564 0.219809 0.238974 0.212109 0.203512 0.262628 0.139122 0.209955 0.25898 0.121253 0.146258 0.252119 0.189315 0.348857 0.192001 0.177253 0.18511 0.318837 0.17274 0.241215 0.234717 0.151474 0.197845 0.172184 0.177052 0.188133 0.181487 0.17894 0.391196 0.239734 0.201666 0.155442 0.284127 0.361163 0.176645 0.213435 0.144676 0.178563 0.184529 0.194893 0.16374 0.194494 0.248919 0.199966 0.294091 0.192567 0.137726 0.205497 0.377547 0.160435 0.177296 0.270629 0.198489 0.277592 0.225705 0.181984 0.208872 0.153021 0.217286 0.250653 0.232345 0.227382 0.207033 0.15112 0.171822 0.200652 0.230417 0.190486 0.226978 0.221965 0.265418 0.191157 0.265547 0.158963 0.233287 0.2337 0.208065 0.154535 0.264492 0.21634 0.242261 0.179481 0.146798 0.319285 0.173971 0.182623 0.183595 0.256185 0.164849 0.237656 0.215889 0.191033 0.238991 0.184443 0.195366 0.305512 0.235429 0.228127 0.236541 0.245111 0.229651 0.230195 0.155211 0.241682 0.216595 0.163756 0.210027 0.170744 0.208308 0.190324 0.159395 0.2384 0.215245 0.264588 0.167321 0.311468 0.228536 0.189227 0.21133 0.201602 0.220131 0.225132 0.186355 0.23987 0.373877 0.169155 0.270835 0.255913 0.191739 0.243021 0.191792 0.25374 0.171217 0.159751 0.202927 0.309625 0.177031 0.206689 0.216633 0.288582 0.282464 0.187681 0.18189 0.229226 0.187604 0.25857 0.275821 0.194373 0.147093 0.214847 0.149519 0.205432 0.203453 0.221115 0.180499 0.253979 0.31916 0.165059 0.221428 0.137058 0.221507 0.26111 0.213237 0.327934 0.139098 0.16988 0.212053 0.273156 0.223122 0.212503 0.213628 0.20651 0.220121 0.189289 0.134708 0.21263 0.337963 0.233443 0.180724 0.250797 0.191758 0.151158 0.193584 0.187774 0.160854 0.19613 0.387155 0.331756 0.227498 0.210868 0.238764 0.208559 0.174363 0.195667 0.250454 0.187979 0.214151 0.22076 0.152997 0.185985 0.187072 0.246929 0.422118 0.106415 0.344988 0.116554 0.252099 0.203194 0.22397 0.349675 0.169378 0.354585 0.230294 0.173531 0.199079 0.232905 0.208504 0.165396 0.252412 0.183778 0.193563 0.241733 0.208057 0.178622 0.222577 0.182392 0.267264 0.259241 0.181699 0.227966 0.240612 0.197367 0.189627 0.128022 0.207789 0.229405 0.227657 0.392709 0.204927 0.188977 0.174428 0.162605 0.198213 0.20212 0.311971 0.258325 0.189808 0.168225 0.211795 0.302341 0.156814 0.218536 0.185633 0.239393 0.245809 0.280321 0.227455 0.24768 0.159827 0.264748 0.197062 0.166835 0.1579 0.161183 0.218967 0.226432 0.156546 0.231311 0.176592 0.262075 0.18962 0.201337 0.307211 0.208688 0.231602 0.359707 0.16772 0.245945 0.327982 0.218915 0.172997 0.188598 0.164793 0.212792 0.236814 0.189839 0.170886 0.212633 0.324881 0.218563 0.209498 0.204217 0.195639 0.247268 0.187592 0.168787 0.171639 0.150795 0.213918 0.382474 0.283452 0.213796 0.233259 0.1739 0.199426 0.224043 0.446501 0.154955 0.190573 0.25372 0.170737 0.1899 0.187603 0.197291 0.251914 0.177571 0.16152 0.214652 0.22332 0.187381 0.257945 0.376001 0.188233 0.241144 0.225396 0.19618 0.147204 0.187352 0.202311 0.24281 0.197581 0.31419 0.250592 0.115361 0.263892 0.155547 0.180945 0.222929 0.1705 0.230904 0.203517 0.227021 0.16291 0.204403 0.195474 0.268045 0.226956 0.16908 0.151408 0.202749 0.19765 0.242888 0.172872 0.192409 0.219095 0.31887 0.152231 0.186755 0.23601 0.215834 0.233374 0.160131 0.204254 0.278227 0.280624 0.208589 0.192215 0.167925 0.156208 0.226438 0.333329 0.205312 0.394569 0.499543 0.225985 0.197851 0.224949 0.231311 0.264459 0.203908 0.233545 0.16896 0.348186 0.206507 0.220755 0.160865 0.1944 0.227322 0.440547 0.351366 0.185278 0.185183 0.22704 0.217496 0.208261 0.255317 0.255393 0.13817 0.162529 0.266165 0.195381 0.243408 0.221481 0.230747 0.179377 0.168538 0.219833 0.181793 0.170123 0.175698 0.153836 0.278811 0.374566 0.186939 0.177772 0.214835 0.196381 0.228222 0.250217 0.279282 0.220819 0.227097 0.270474 0.183295 0.250454 0.205413 0.175691 0.312591 0.221383 0.395739 0.312435 0.170681 0.197736 0.168583 0.197877 0.291335 0.17791 0.252802 0.203421 0.221956 0.277748 0.218761 0.192907 0.262229 0.195982 0.284449 0.145633 0.243801 0.24442 0.221999 0.312926 0.241497 0.201007 0.302999 0.21086 0.203566 0.277272 0.181857 0.172218 0.236629 0.182372 0.284157 0.272622 0.225269 0.197008 0.250537 0.210721 0.17952 0.270447 0.220893 0.258048 0.158512 0.231829 0.272997 0.157319 0.168473 0.244472 0.225119 0.202538 0.25498 0.214426 0.205963 0.199192 0.192237 0.211226 0.196844 0.187147 0.167375 0.204043 0.310633 0.219111 0.20002 0.179792 0.209861 0.222578 0.179712 0.161976 0.204327 0.310482 0.202696 0.173593 0.202557 0.154968 0.211989 0.253052 0.312248 0.235669 0.228012 0.170319 0.226037 0.155823 0.323794 0.330386 0.242247 0.197057 0.198924 0.171434 0.23472 0.201159 0.172954 0.335307 0.217496 0.181024 0.219099 0.167395 0.205796 0.2714 0.257132 0.217506 0.317854 0.331396 0.25783 0.190466 0.213451 0.152535 0.209141 0.170084 0.168145 0.220901 0.216539 0.213948 0.161722 0.137255 0.51445 0.168878 0.164768 0.220052 0.218041 0.209833 0.252839 0.193558 0.181193 0.190316 0.245163 0.259084 0.236704 0.194804 0.20751 0.143671 0.178263 0.17909 0.196044 0.32529 0.251026 0.277489 0.233825 0.269091 0.223747 0.279871 0.21079 0.217189 0.349653 0.19186 0.228004 0.211953 0.178762 0.23843 0.216873 0.23964 0.189913 0.264208 0.251239 0.156337 0.22064 0.174541 0.226594 0.336344 0.225812 0.189841 0.152799 0.279935 0.216876 0.229113 0.200344 0.168261 0.1515 0.180628 0.234959 0.227595 0.180516 0.210599 0.180549 0.184988 0.20505 0.205242 0.20802 0.184056 0.17997 0.247108 0.191051 0.224323 0.505734 0.183319 0.224556 0.170376 0.445638 0.241509 0.206027 0.304087 0.246093 0.186353 0.219425 0.181789 0.18959 0.238821 0.229441 0.215655 0.233767 0.338171 0.208629 0.20028 0.215128 0.307776 0.196185 0.167519 0.220951 0.188588 0.187153 0.184394 0.174089 0.207032 0.243532 0.259983 0.237466 0.196898 0.194766 0.201618 0.28287 0.251813 0.224254 0.166095 0.282685 0.199957 0.218241 0.187594 0.296093 0.210137 0.20762 0.270575 0.26017 0.211321 0.177226 0.149355 0.258134 0.21814 0.201425 0.293922 0.568611 0.295506 0.211349 0.170565 0.198215 0.211363 0.279449 0.310756 0.190406 0.179982 0.29526 0.145965 0.198665 0.17355 0.225618 0.198843 0.257595 0.218682 0.219216 0.168171 0.245022 0.279247 0.254746 0.307589 0.187613 0.226769 0.242956 0.206057 0.159609 0.184053 0.143681 0.249347 0.22425 0.188352 0.235941 0.194122 0.26051 0.210049 0.1817 0.220168 0.202591 0.257684 0.160623 0.376192 0.170598 0.165106 0.15927 0.186598 0.278064 0.262584 0.165841 0.223583 0.183595 0.161967 0.222596 0.164135 0.21108 0.271042 0.212486 0.175381 0.21602 0.165677 0.263014 0.214321 0.181313 0.217124 0.219385 0.262021 0.202169 0.234638 0.189105 0.189134 0.331813 0.168405 0.156222 0.322208 0.215112 0.183436 0.177441 0.146784 0.303849 0.216556 0.191912 0.204518 0.243673 0.329784 0.181442 0.170469 0.291847 0.143284 0.16753 0.140944 0.273807 0.211148 0.169481 0.171658 0.182457 0.099247 0.344651 0.174891 0.233624 0.193752 0.321167 0.209178 0.213388 0.20544 0.173054 0.224743 0.187373 0.220627 0.252146 0.18297 0.194935 0.184989 0.137204 0.17222 0.216998 0.168366 0.204216 0.195239 0.181552 0.181439 0.187082 0.153952 0.32869 0.205574 0.216988 0.235912 0.189647 0.233765 0.13935 0.212406 0.340976 0.16274 0.282618 0.215032 0.253185 0.1571 0.194609 0.277801 0.185466 0.174802 0.279586 0.192481 0.202636 0.199436 0.267234 0.268719 0.188847 0.307071 0.235268 0.217078 0.191929 0.150674 0.225292 0.177312 0.159417 0.220078 0.208678 0.16619 0.167587 0.219145 0.162184 0.18755 0.197153 0.200582 0.341263 0.140527 0.216683 0.173795 0.268461 0.163528 0.37339 0.301437 0.186144 0.200729 0.156643 0.253806 0.208154 0.234934 0.202839 0.290611 0.232148 0.250699 0.276606 0.29556 0.190251 0.248114 0.147264 0.216552 0.198107 0.239066 0.228653 0.190881 0.243586 0.138756 0.213119 0.235716 0.19785 0.113691 0.176615 0.195127 0.224503 0.27389 0.192991 0.16982 0.301443 0.243976 0.171245 0.210483 0.206787 0.230388 0.2292 0.147426 0.183411 0.180082 0.211088 0.164719 0.208282 0.234761 0.190826 0.186094 0.362261 0.130066 0.176169 0.213626 0.307393 0.181525 0.231458 0.231428 0.28858 0.120082 0.278366 0.181593 0.209308 0.239985 0.219727 0.23101 0.258567 0.151317 0.176744 0.157823 0.17485 0.196612 0.206978 0.235078 0.276803 0.250286 0.199888 0.159089 0.197719 0.186858 0.196891 0.230503 0.228479 0.194937 0.143275 0.182655 0.254217 0.239347 0.263629 0.245069 0.203328 0.219396 0.242226 0.332494 0.307093 0.214993 0.148866 0.246008 0.289999 0.161221 0.262341 0.27531 0.169127 0.376075 0.199306 0.237387 0.137877 0.205165 0.408509 0.173003 0.193972 0.208691 0.198339 0.20013 0.170464 0.278355 0.226973 0.182068 0.198523 0.242829 0.212875 0.335446 0.277698 0.215893 0.276681 0.181049 0.208186 0.242637 0.213253 0.213592 0.173112 0.25074 0.192978 0.230149 0.269885 0.259714 0.153693 0.156405 0.183263 0.24421 0.170866 0.192807 0.316486 0.260069 0.19414 0.129796 0.192197 0.253264 0.149589 0.229527 0.271938 0.247946 0.365825 0.178131 0.218047 0.174941 0.302977 0.208841 0.216272 0.211365 0.227409 0.141777 0.190281 0.156575 0.179368 0.244756 0.240108 0.370096 0.188118 0.220804 0.230333 0.247989 0.26134 0.17729 0.198413 0.178397 0.211252 0.183739 0.201152 0.205616 0.323547 0.198617 0.150904 0.171332 0.268078 0.175102 0.238083 0.208334 0.195694 0.280975 0.230256 0.242902 0.197819 0.204517 0.194574 0.190731 0.36237 0.197161 0.19449 0.180017 0.196078 0.224389 0.202459 0.234853 0.238181 0.210833 0.240359 0.229272 0.216201 0.152928 0.218109 0.277874 0.346381 0.194394 0.18626 0.188062 0.327204 0.144859 0.232843 0.171366 0.183936 0.190749 0.177405 0.22606 0.172744 0.23683 0.301163 0.414293 0.211337 0.219869 0.211175 0.184881 0.184771 0.215937 0.246162 0.297166 0.194441 0.282912 0.186541 0.20433 0.190816 0.24565 0.293597 0.197892 0.176691 0.211929 0.234598 0.26874 0.188825 0.289382 0.180266 0.297202 0.161616 0.200941 0.191616 0.170613 0.195972 0.143642 0.256814 0.125352 0.180669 0.195445 0.230906 0.172723 0.221086 0.185193 0.204739 0.241652 0.214425 0.330573 0.364922 0.20525 0.276545 0.165466 0.36864 0.176448 0.245331 0.213111 0.200123 0.262577 0.171518 0.23376 0.119244 0.280563 0.198711 0.222224 0.231828 0.301398 0.211656 0.162635 0.177042 0.340298 0.229334 0.216574 0.345853 0.148967 0.197498 0.194825 0.184554 0.151946 0.17857 0.185107 0.19653 0.235294 0.216762 0.218521 0.12383 0.279246 0.251441 0.210006 0.146587 0.128376 0.185576 0.225939 0.159182 0.320646 0.143072 0.194932 0.262782 0.228906 0.192957 0.183348 0.248739 0.175389 0.211233 0.238279 0.282187 0.234628 0.1467 0.294617 0.326503 0.134931 0.205552 0.194762 0.312941 0.158796 0.183553 0.184759 0.248451 0.183784 0.212856 0.224554 0.218423 0.22578 0.210339 0.203489 0.287066 0.188176 0.291427 0.17686 0.169038 0.175795 0.177522 0.195336 0.233433 0.21489 0.151712 0.270438 0.213113 0.143738 0.16366 0.184074 0.156223 0.17686 0.276277 0.219923 0.211711 0.257255 0.275854 0.207692 0.388379 0.225055 0.185971 0.173915 0.203891 0.398087 0.192057 0.172525 0.251942 0.232589 0.260697 0.219992 0.20522 0.165365 0.276189 0.221422 0.216057 0.192348 0.202856 0.200787 0.181968 0.15576 0.224888 0.254018 0.248581 0.273766 0.28598 0.218597 0.279891 0.265559 0.227598 0.213856 0.201635 0.247724 0.180687 0.195987 0.19506 0.186911 0.207388 0.188235 0.301338 0.1829 0.25865 0.219215 0.20038 0.379174 0.198097 0.170395 0.201044 0.20342 0.177506 0.18128 0.223003 0.240808 0.217703 0.211127 0.17936 0.181563 0.29615 0.195203 0.189915 0.240809 0.23717 0.212809 0.19644 0.173598 0.247632 0.215067 0.231184 0.146446 0.248358 0.159744 0.240479 0.266187 0.198386 0.187093 0.253738 0.187371 0.246746 0.185911 0.176265 0.259215 0.21357 0.195424 0.189828 0.16076 0.190111 0.196158 0.348731 0.425219 0.219865 0.223061 0.234964 0.190476 0.293156 0.284626 0.208864 0.220283 0.238751 0.198906 0.215444 0.152216 0.186431 0.183113 0.19723 0.163489 0.243372 0.190056 0.240335 0.153688 0.264403 0.118294 0.188796 0.432532 0.267456 0.36688 0.221817 0.220856 0.202154 0.219521 0.162097 0.222316 0.225669 0.247323 0.190754 0.354435 0.342488 0.215826 0.375593 0.216714 0.234711 0.298743 0.196714 0.237568 0.270703 0.240921 0.209058 0.20709 0.185641 0.215538 0.299779 0.182044 0.277206 0.171484 0.256935 0.235834 0.236484 0.151127 0.259121 0.148356 0.20209 0.236712 0.210217 0.158046 0.246565 0.314336 0.266594 0.142491 0.244733 0.191507 0.185805 0.291143 0.175395 0.183058 0.191714 0.128669 0.196771 0.186208 0.250229 0.224718 0.201347 0.151723 0.29284 0.138252 0.211024 0.314624 0.203801 0.213051 0.252377 0.209594 0.182219 0.133487 0.19404 0.39339 0.197218 0.174167 0.333195 0.24461 0.184624 0.248033 0.224357 0.253577 0.340173 0.161713 0.194246 0.223892 0.162263 0.206744 0.231191 0.166203 0.228167 0.195092 0.1973 0.200967 0.127889 0.216436 0.280223 0.19269 0.146241 0.126488 0.182065 0.226024 0.228594 0.210324 0.375505 0.180998 0.245256 0.188775 0.193472 0.185426 0.213729 0.202388 0.185696 0.313863 0.281444 0.210197 0.181662 0.225365 0.240246 0.212856 0.212039 0.279764 0.176979 0.190357 0.244739 0.191608 0.226983 0.206895 0.168401 0.164831 0.199246 0.381433 0.18331 0.314064 0.219133 0.166521 0.212745 0.221101 0.20643 0.228425 0.166781 0.224454 0.231818 0.188087 0.219348 0.208976 0.323424 0.230664 0.341256 0.221316 0.175988 0.281625 0.153294 0.263242 0.223245 0.193771 0.312292 0.239469 0.158025 0.181813 0.178154 0.154447 0.181416 0.233754 0.225637 0.200649 0.150712 0.171892 0.224952 0.202711 0.170137 0.171836 0.169325 0.234668 0.213953 0.218277 0.229906 0.208748 0.180643 0.259106 0.187427 0.211417 0.19224 0.243133 0.24971 0.302524 0.214354 0.233155 0.208339 0.190582 0.224484 0.213727 0.253671 0.182313 0.198217 0.206605 0.250053 0.180193 0.240035 0.166109 0.188659 0.192821 0.190994 0.168596 0.164122 0.202237 0.238072 0.145363 0.250562 0.212073 0.264763 0.267798 0.184631 0.274526 0.143698 0.185622 0.240448 0.187405 0.230093 0.216945 0.2768 0.223438 0.160215 0.662305 0.209937 0.159769 0.186272 0.298257 0.311477 0.166087 0.208241 0.206002 0.240633 0.211116 0.373898 0.220763 0.220603 0.211437 0.222755 0.210373 0.278439 0.230647 0.203798 0.173399 0.197339 0.223866 0.179816 0.176918 0.160345 0.203568 0.275943 0.190057 0.200562 0.160133 0.208376 0.20651 0.276842 0.134329 0.1641 0.194955 0.3966 0.214825 0.293561 0.221986 0.216033 0.284414 0.200415 0.23446 0.315146 0.221045 0.422361 0.171654 0.195081 0.295145 0.193501 0.245567 0.260379 0.177816 0.23766 0.20067 0.208197 0.198916 0.227357 0.188488 0.344142 0.223653 0.201177 0.177033 0.143487 0.219831 0.134802 0.1542 0.204522 0.148403 0.164446 0.24756 0.2214 0.202534 0.182203 0.186749 0.210957 0.156229 0.225164 0.220126 0.181926 0.359802 0.204581 0.188135 0.193951 0.217961 0.207434 0.184312 0.176296 0.187511 0.22256 0.210407 0.153888 0.212121 0.190297 0.17193 0.370311 0.183742 0.249608 0.196264 0.153818 0.190617 0.162675 0.200461 0.250517 0.160508 0.216874 0.132658 0.340792 0.230344 0.255154 0.278189 0.234787 0.266288 0.220006 0.215879 0.172501 0.194912 0.211343 0.21928 0.192485 0.190975 0.232451 0.19422 0.214124 0.148634 0.253469 0.260599 0.186228 0.30533 0.151699 0.260337 0.202372 0.205418 0.18147 0.191342 0.215484 0.219941 0.207516 0.227725 0.301892 0.199602 0.158271 0.188914 0.258198 0.275153 0.245872 0.230557 0.266246 0.199454 0.239516 0.387285 0.271389 0.29542 0.223406 0.240367 0.173695 0.239426 0.209042 0.183823 0.177195 0.236083 0.233903 0.223958 0.159939 0.199829 0.175668 0.201148 0.212649 0.183532 0.16669 0.297363 0.347608 0.241704 0.158502 0.23815 0.238526 0.241765 0.161428 0.223081 0.189774 0.210065 0.194456 0.170777 0.154102 0.224568 0.237687 0.198255 0.244054 0.324425 0.200145 0.266295 0.196388 0.206405 0.194665 0.267003 0.180725 0.243871 0.360161 0.223389 0.147078 0.282675 0.227368 0.147247 0.248641 0.25552 0.264928 0.20016 0.386816 0.198218 0.2327 0.233233 0.194062 0.332074 0.245238 0.233491 0.244208 0.179534 0.22129 0.238342 0.176331 0.177872 0.287729 0.291555 0.170451 0.303618 0.25496 0.213006 0.215105 0.301192 0.293717 0.163829 0.215045 0.321837 0.302417 0.296647 0.204793 0.121307 0.278043 0.210617 0.298234 0.262611 0.199216 0.237051 0.196276 0.206126 0.180423 0.178651 0.273833 0.207939 0.142614 0.199641 0.272817 0.249967 0.278943 0.201564 0.286156 0.152845 0.189684 0.189668 0.234045 0.162305 0.205949 0.163487 0.278869 0.136652 0.288122 0.173614 0.155428 0.195908 0.154735 0.261974 0.198461 0.220525 0.194289 0.252276 0.158885 0.181897 0.173468 0.179473 0.248717 0.374407 0.229399 0.203235 0.196816 0.304573 0.225977 0.252386 0.209445 0.153218 0.204556 0.241856 0.191768 0.227338 0.272526 0.23602 0.222297 0.178539 0.191017 0.237715 0.205033 0.146854 0.175215 0.245157 0.226018 0.189426 0.202047 0.183544 0.273209 0.216036 0.21232 0.233591 0.172066 0.25208 0.342007 0.12502 0.261658 0.231996 0.229274 0.170907 0.194555 0.30377 0.214558 0.208761 0.177524 0.193478 0.205597 0.214887 0.168678 0.236767 0.204643 0.175617 0.375601 0.110495 0.221253 0.21662 0.235651 0.141697 0.158069 0.204275 0.197539 0.150483 0.226127 0.193966 0.236632 0.164065 0.185847 0.309858 0.172838 0.187991 0.160358 0.103359 0.144089 0.180415 0.177997 0.277542 0.226563 0.225151 0.173064 0.209619 0.240054 0.163193 0.184265 0.176755 0.230768 0.21646 0.187716 0.260321 0.210238 0.205703 0.247339 0.181933 0.126445 0.199416 0.15515 0.221554 0.200302 0.223929 0.172531 0.218512 0.256329 0.202 0.130368 0.139097 0.23867 0.264748 0.242223 0.410433 0.413099 0.164701 0.292658 0.222256 0.258805 0.196617 0.229698 0.186819 0.178395 0.303453 0.189295 0.21716 0.254604 0.175328 0.405073 0.220965 0.195518 0.225033 0.176593 0.406685 0.197782 0.170986 0.206681 0.210123 0.151272 0.189696 0.216532 0.185098 0.115726 0.087921 0.162177 0.200464 0.204637 0.213001 0.141352 0.193551 0.197534 0.302396 0.181888 0.218075 0.158762 0.161016 0.207227 0.205 0.171758 0.236699 0.193441 0.210142 0.206992 0.172827 0.196477 0.216707 0.161231 0.184225 0.171599 0.255937 0.265495 0.288286 0.227324 0.251383 0.233874 0.198781 0.200179 0.201781 0.172971 0.285973 0.305448 0.306708 0.257573 0.290808 0.216293 0.125289 0.298235 0.225293 0.219173 0.221756 0.219794 0.204746 0.185196 0.190752 0.216971 0.140851 0.141079 0.201121 0.166175 0.275524 0.207654 0.24481 0.249109 0.209646 0.158995 0.179333 0.264283 0.264825 0.148654 0.168816 0.205431 0.195362 0.206728 0.265333 0.190383 0.17255 0.265402 0.305415 0.142861 0.177852 0.155254 0.256299 0.121774 0.216692 0.19146 0.178792 0.286078 0.225074 0.2754 0.226635 0.251459 0.29153 0.186526 0.182722 0.316759 0.194549 0.324015 0.196706 0.229819 0.152482 0.279673 0.190269 0.203285 0.218078 0.266586 0.212927 0.176571 0.149923 0.186065 0.281359 0.169434 0.221411 0.198838 0.291173 0.262003 0.175323 0.277022 0.248002 0.220206 0.207556 0.229182 0.183638 0.231776 0.326672 0.221761 0.221243 0.22649 0.172126 0.241313 0.334846 0.199587 0.201097 0.185716 0.207681 0.192715 0.216822 0.232901 0.163732 0.203534 0.253921 0.204371 0.207817 0.164905 0.133038 0.213946 0.200384 0.191607 0.20374 0.388166 0.179256 0.261552 0.390182 0.242797 0.193508 0.381035 0.215218 0.2375 0.178937 0.199593 0.170413 0.256926 0.250883 0.149438 0.189503 0.16145 0.240987 0.215875 0.191241 0.229156 0.225381 0.200593 0.217476 0.162239 0.226808 0.192402 0.232375 0.361678 0.16424 0.200737 0.205996 0.159645 0.324869 0.348212 0.225733 0.223032 0.253458 0.167454 0.228957 0.192601 0.2179 0.227955 0.280298 0.148509 0.20131 0.26152 0.263747 0.225687 0.167408 0.137769 0.203269 0.15862 0.332561 0.215785 0.243535 0.282886 0.186909 0.190277 0.247122 0.138787 0.228069 0.191969 0.248288 0.155765 0.201241 0.177045 0.231088 0.218787 0.201321 0.26259 0.189927 0.247755 0.329475 0.296966 0.224381 0.279427 0.216288 0.22353 0.21967 0.204164 0.213464 0.207568 0.175772 0.185519 0.189597 0.183359 0.196484 0.207552 0.168072 0.208961 0.208749 0.165986 0.264267 0.230792 0.246247 0.310964 0.203983 0.190677 0.287072 0.21364 0.217056 0.227899 0.132911 0.269793 0.228903 0.14729 0.241343 0.165992 0.271494 0.222889 0.208995 0.212914 0.215072 0.176528 0.193178 0.195058 0.181732 0.239693 0.271194 0.10724 0.141055 0.223415 0.202582 0.249853 0.2755 0.231924 0.196459 0.278029 0.201916 0.168638 0.202407 0.232381 0.249473 0.214373 0.211552 0.242326 0.22353 0.317031 0.155974 0.144757 0.234652 0.136254 0.192208 0.203084 0.181575 0.145721 0.269559 0.267964 0.165433 0.111417 0.286419 0.191068 0.173711 0.215723 0.171494 0.212627 0.276836 0.399073 0.256116 0.293996 0.208954 0.215933 0.20973 0.171288 0.259672 0.189076 0.19374 0.182657 0.143422 0.12519 0.194295 0.215 0.29657 0.170177 0.352467 0.197387 0.235542 0.19138 0.223155 0.181389 0.197095 0.132406 0.157641 0.180485 0.215874 0.180045 0.166473 0.306402 0.153394 0.276226 0.3755 0.216072 0.248882 0.206489 0.141162 0.142962 0.222038 0.200139 0.192318 0.151331 0.238679 0.215196 0.191446 0.225417 0.282933 0.141595 0.189248 0.198611 0.141576 0.271069 0.191289 0.210268 0.195444 0.164127 0.203022 0.25583 0.212089 0.1714 0.271155 0.299391 0.250588 0.184843 0.224795 0.412182 0.210704 0.185782 0.170203 0.261779 0.223184 0.157015 0.167582 0.199859 0.096694 0.184265 0.215447 0.143406 0.266372 0.143217 0.190264 0.247006 0.149029 0.24052 0.276259 0.133227 0.265557 0.208212 0.227624 0.213058 0.185862 0.219328 0.213341 0.202572 0.275436 0.172598 0.149075 0.253344 0.255814 0.128999 0.252551 0.190066 0.187986 0.156015 0.221614 0.189283 0.331532 0.17793 0.139829 0.268387 0.351852 0.16467 0.178319 0.16857 0.288686 0.317991 0.147125 0.161078 0.203963 0.169944 0.245208 0.278023 0.169014 0.17326 0.138978 0.247387 0.238095 0.26181 0.244028 0.251919 0.172195 0.14672 0.291942 0.258393 0.20011 0.145547 0.320137 0.170686 0.254234 0.240599 0.252496 0.210686 0.202613 0.20574 0.203701 0.182476 0.203502 0.197242 0.171826 0.31041 0.22642 0.182345 0.187584 0.160911 0.145584 0.102482 0.187043 0.361409 0.217509 0.220119 0.256818 0.183619 0.163328 0.190966 0.257568 0.151202 0.23619 0.217455 0.257581 0.181586 0.206458 0.327925 0.263439 0.267406 0.148313 0.233107 0.168432 0.169475 0.233779 0.200847 0.146452 0.233941 0.255763 0.17434 0.237573 0.189206 0.189257 0.292318 0.230353 0.273628 0.214819 0.229092 0.159989 0.283012 0.267308 0.201899 0.382227 0.257008 0.255658 0.248371 0.180941 0.193341 0.16972 0.231973 0.263314 0.23578 0.171409 0.191223 0.206732 0.234335 0.206254 0.204844 0.188823 0.206864 0.217859 0.243186 0.209768 0.277135 0.216274 0.299705 0.14824 0.18806 0.265025 0.186799 0.232862 0.17256 0.14803 0.228107 0.220481 0.24384 0.180517 0.312066 0.182633 0.179889 0.208609 0.234989 0.217033 0.187809 0.12922 0.217165 0.316146 0.209963 0.161255 0.174565 0.234458 0.331032 0.218013 0.193879 0.177369 0.22104 0.216568 0.195072 0.268008 0.232719 0.242865 0.198699 0.248178 0.142592 0.183655 0.224911 0.210949 0.236958 0.132951 0.177617 0.208369 0.228276 0.196572 0.220355 0.167859 0.249046 0.227484 0.264682 0.343287 0.201813 0.233824 0.266497 0.278483 0.179299 0.193622 0.230032 0.349756 0.241726 0.257899 0.149428 0.20639 0.255397 0.268343 0.220671 0.180977 0.181267 0.244474 0.185073 0.199097 0.257047 0.238903 0.201511 0.211199 0.259529 0.226913 0.201008 0.22094 0.267348 0.22183 0.19212 0.199965 0.259447 0.210665 0.148992 0.219715 0.236723 0.200841 0.251319 0.167928 0.180745 0.192313 0.349945 0.194596 0.17102 0.205728 0.195307 0.157981 0.2636 0.280611 0.263757 0.166578 0.237379 0.17909 0.212031 0.204296 0.220863 0.307974 0.191868 0.345015 0.170636 0.188675 0.292823 0.393025 0.238731 0.299044 0.175712 0.191075 0.231464 0.272837 0.198575 0.230535 0.200211 0.271289 0.205748 0.202026 0.192523 0.207093 0.203182 0.219654 0.214837 0.214511 0.284836 0.203995 0.155487 0.262495 0.159975 0.218996 0.221193 0.20361 0.229606 0.1932 0.218153 0.195442 0.200592 0.344549 0.211521 0.314196 0.221889 0.260979 0.185552 0.149416 0.196817 0.185683 0.22306 0.149353 0.182266 0.210061 0.228337 0.208321 0.232709 0.20833 0.196575 0.158122 0.160798 0.198104 0.279958 0.272036 0.181481 0.174543 0.22555 0.1645 0.224203 0.207017 0.149355 0.237567 0.185914 0.162651 0.19543 0.198954 0.183701 0.201869 0.219999 0.185485 0.173206 0.127757 0.245032 0.205662 0.206605 0.196818 0.23228 0.223298 0.229162 0.158992 0.232256 0.243296 0.244211 0.212443 0.218799 0.174939 0.187116 0.207696 0.260837 0.22316 0.161933 0.218968 0.209022 0.228571 0.195959 0.184866 0.181556 0.214136 0.136115 0.218041 0.236878 0.220319 0.496668 0.185957 0.249973 0.299266 0.213082 0.239479 0.192375 0.217301 0.191042 0.305693 0.281449 0.282444 0.308368 0.183742 0.210923 0.290723 0.18014 0.196221 0.180897 0.209403 0.273449 0.192027 0.201848 0.179589 0.247367 0.194427 0.292273 0.247951 0.207063 0.218915 0.215864 0.217133 0.2295 0.242278 0.197245 0.272272 0.222072 0.355124 0.206777 0.249415 0.260962 0.219599 0.250576 0.16671 0.300081 0.244223 0.232322 0.190382 0.195455 0.191103 0.137093 0.277923 0.120976 0.203597 0.252089 0.231045 0.16835 0.44148 0.237384 0.242242 0.271289 0.205709 0.435127 0.186267 0.23489 0.245281 0.160841 0.148974 0.27719 0.260334 0.219419 0.206522 0.180556 0.246507 0.200315 0.217533 0.19684 0.323496 0.265772 0.189359 0.156651 0.33758 0.254662 0.179702 0.249596 0.200663 0.313212 0.133896 0.237429 0.284399 0.194866 0.21511 0.228735 0.243159 0.212144 0.299901 0.155514 0.189949 0.188491 0.158634 0.193973 0.25714 0.250529 0.178708 0.186705 0.151247 0.295488 0.229334 0.176464 0.201301 0.231715 0.185786 0.244011 0.142876 0.226238 0.254335 0.312798 0.335175 0.206609 0.197597 0.147209 0.217768 0.270163 0.243256 0.224545 0.214693 0.252148 0.201783 0.210279 0.216933 0.219194 0.217186 0.169531 0.135378 0.191888 0.214571 0.182229 0.257357 0.200464 0.266585 0.129532 0.173526 0.204687 0.230292 0.225179 0.198206 0.215994 0.149489 0.45024 0.222131 0.330527 0.218053 0.224737 0.291447 0.390478 0.224373 0.25474 0.221151 0.2068 0.197066 0.216378 0.208701 0.224441 0.229879 0.235203 0.208687 0.206128 0.204199 0.200374 0.224048 0.428724 0.205479 0.251155 0.17513 0.20795 0.228378 0.203219 0.196233 0.210169 0.282007 0.212127 0.272495 0.254305 0.160637 0.148173 0.227633 0.221321 0.165865 0.197069 0.197712 0.239311 0.182537 0.205405 0.222615 0.1954 0.295618 0.196271 0.147267 0.229285 0.205741 0.334702 0.197723 0.2753 0.205849 0.229976 0.225517 0.185494 0.279212 0.22357 0.17807 0.21247 0.199477 0.207294 0.239979 0.230406 0.259167 0.163803 0.207054 0.313422 0.147521 0.208605 0.179376 0.308488 0.1953 0.182029 0.183843 0.190126 0.1905 0.189922 0.206424 0.220788 0.165396 0.203209 0.264151 0.199109 0.181218 0.260126 0.199883 0.166126 0.333166 0.320918 0.240421 0.212112 0.184157 0.227409 0.236187 0.176978 0.206623 0.253012 0.21193 0.186149 0.205043 0.200782 0.17387 0.175145 0.266527 0.20166 0.435567 0.189009 0.208749 0.173595 0.172979 0.186162 0.380855 0.236104 0.161269 0.173986 0.182237 0.173237 0.28048 0.203659 0.220842 0.186742 0.410841 0.245108 0.137782 0.190508 0.142294 0.331046 0.164042 0.214739 0.17505 0.204327 0.225963 0.161165 0.322968 0.256832 0.133683 0.259873 0.245581 0.188764 0.210478 0.163843 0.209977 0.191253 0.236476 0.421002 0.193476 0.178112 0.337288 0.255152 0.241998 0.196116 0.402636 0.182578 0.359509 0.233217 0.21401 0.149955 0.202827 0.241121 0.225551 0.218476 0.201498 0.245009 0.380938 0.222873 0.182126 0.185101 0.253995 0.335816 0.216935 0.183905 0.239519 0.231783 0.210684 0.236893 0.189911 0.211222 0.227075 0.237142 0.224551 0.219545 0.370251 0.214133 0.18479 0.152638 0.278886 0.202336 0.169557 0.172429 0.207358 0.210304 0.214361 0.189416 0.235263 0.248659 0.223007 0.266902 0.224855 0.217856 0.222367 0.154007 0.169462 0.164217 0.184057 0.189316 0.207798 0.30467 0.207789 0.183414 0.170808 0.180536 0.181972 0.172912 0.240683 0.268181 0.266238 0.244039 0.204921 0.169059 0.152624 0.255604 0.211522 0.283369 0.195872 0.211102 0.201849 0.236597 0.222753 0.168337 0.181861 0.308274 0.179387 0.287574 0.178725 0.309162 0.229476 0.208125 0.220366 0.132126 0.229946 0.201584 0.255485 0.255372 0.275834 0.210338 0.223143 0.211101 0.187387 0.224879 0.207506 0.269005 0.253684 0.160644 0.228056 0.386271 0.287701 0.265768 0.197482 0.21786 0.316047 0.122206 0.157482 0.218939 0.208634 0.192866 0.245738 0.17862 0.203004 0.293715 0.211432 0.195005 0.212218 0.176508 0.198637 0.228432 0.274271 0.178932 0.206068 0.207581 0.201263 0.259716 0.226635 0.228264 0.270557 0.170966 0.232614 0.197397 0.221364 0.212765 0.209078 0.193293 0.186895 0.209061 0.205549 0.224936 0.205191 0.220635 0.176712 0.161075 0.269506 0.247706 0.193656 0.231065 0.204924 0.097186 0.496358 0.182975 0.286954 0.180575 0.166999 0.141959 0.236955 0.178438 0.144333 0.17723 0.176312 0.190184 0.228494 0.173073 0.145733 0.209474 0.140284 0.180644 0.170297 0.32745 0.274855 0.164161 0.190131 0.195796 0.237804 0.227271 0.201394 0.257129 0.210625 0.242051 0.207623 0.208909 0.213069 0.166074 0.183847 0.163358 0.147264 0.208877 0.206973 0.147075 0.196396 0.170125 0.127129 0.154176 0.220456 0.193334 0.335892 0.256591 0.220318 0.226337 0.162797 0.259456 0.252413 0.172386 0.235362 0.151323 0.227622 0.138062 0.158075 0.166806 0.174862 0.185799 0.242063 0.180926 0.240663 0.202187 0.239269 0.173738 0.229319 0.304376 0.205283 0.15453 0.251232 0.263624 0.194293 0.133388 0.213053 0.313943 0.136735 0.216777 0.204275 0.172502 0.205186 0.173257 0.237856 0.228985 0.184235 0.225706 0.346598 0.181456 0.201543 0.163241 0.266429 0.309715 0.284654 0.158562 0.19763 0.164653 0.153368 0.225679 0.367141 0.229787 0.194182 0.195372 0.170137 0.158924 0.183077 0.212842 0.35523 0.179698 0.249253 0.284839 0.284106 0.210512 0.12671 0.19093 0.286818 0.2348 0.199474 0.380104 0.179774 0.230039 0.287283 0.166679 0.201471 0.156638 0.211734 0.129828 0.277657 0.313387 0.239496 0.200165 0.248965 0.253718 0.168315 0.180935 0.314152 0.353271 0.188535 0.195821 0.201161 0.271232 0.227073 0.414963 0.154678 0.198436 0.280302 0.18728 0.239201 0.162038 0.232906 0.168201 0.229457 0.18014 0.17698 0.228194 0.210015 0.237163 0.267643 0.267826 0.208242 0.16196 0.222542 0.244308 0.248842 0.146755 0.240873 0.271258 0.190299 0.151343 0.246426 0.265629 0.221877 0.187861 0.186299 0.142723 0.243845 0.244362 0.195748 0.169615 0.169282 0.171318 0.256268 0.302884 0.193838 0.185645 0.242465 0.276129 0.252075 0.21539 0.179695 0.207448 0.202351 0.314429 0.265459 0.173031 0.233685 0.193769 0.162563 0.249657 0.209806 0.269742 0.185239 0.3788 0.245998 0.272187 0.197901 0.223575 0.150716 0.248642 0.235126 0.256453 0.165963 0.231839 0.204677 0.185202 0.190888 0.244632 0.204429 0.130026 0.30682 0.217536 0.186443 0.206357 0.315692 0.369554 0.256529 0.337116 0.254918 0.273119 0.231227 0.193455 0.354968 0.173862 0.214317 0.230814 0.244683 0.228743 0.168433 0.199529 0.185076 0.205989 0.267294 0.203504 0.235209 0.20178 0.191125 0.218791 0.146304 0.289386 0.146626 0.288865 0.21775 0.199089 0.179599 0.187642 0.298487 0.219816 0.25747 0.198275 0.212506 0.157749 0.167198 0.148715 0.241249 0.280806 0.141991 0.196659 0.315079 0.297418 0.298158 0.403712 0.259883 0.185604 0.193489 0.228837 0.232599 0.269466 0.384549 0.322545 0.250909 0.230822 0.228669 0.150299 0.177763 0.202814 0.238655 0.241314 0.21547 0.303286 0.182228 0.197077 0.29728 0.428164 0.236468 0.227486 0.294125 0.274008 0.311886 0.179162 0.216013 0.185728 0.194541 0.236386 0.260177 0.178824 0.208705 0.384603 0.180971 0.162844 0.200121 0.307404 0.192062 0.33952 0.213439 0.283522 0.210259 0.232033 0.207193 0.304147 0.211711 0.275506 0.299201 0.188994 0.141178 0.205338 0.23981 0.195569 0.229351 0.188046 0.208675 0.244962 0.163363 0.215931 0.138909 0.20091 0.241457 0.245787 0.254568 0.235704 0.226549 0.199444 0.202195 0.210094 0.2954 0.151389 0.222582 0.243245 0.189822 0.258534 0.299322 0.199298 0.209919 0.198324 0.243454 0.188964 0.256726 0.291027 0.155031 0.223907 0.253587 0.19939 0.162156 0.234248 0.209903 0.261452 0.16709 0.213863 0.333725 0.269358 0.259302 0.337769 0.161566 0.174727 0.181864 0.205294 0.189246 0.341172 0.277803 0.255678 0.289393 0.241866 0.227916 0.283016 0.27134 0.223163 0.163995 0.23692 0.157212 0.158332 0.135278 0.328804 0.185372 0.222203 0.152071 0.166204 0.255104 0.174668 0.224907 0.224912 0.211105 0.224271 0.226409 0.235491 0.254596 0.234248 0.20521 0.268136 0.147602 0.179951 0.22937 0.220857 0.258224 0.246447 0.185284 0.293573 0.258709 0.207555 0.173788 0.224855 0.193164 0.178839 0.278704 0.23345 0.204083 0.195419 0.147023 0.172777 0.205046 0.18837 0.24819 0.176196 0.218294 0.205787 0.186831 0.220405 0.250436 0.237711 0.206048 0.146097 0.188547 0.155327 0.243418 0.204877 0.231358 0.232651 0.200977 0.174684 0.182461 0.195217 0.188046 0.235064 0.321073 0.34611 0.164563 0.151514 0.157137 0.196337 0.222884 0.233043 0.207735 0.246581 0.207586 0.19294 0.200726 0.171672 0.148248 0.237108 0.353139 0.216143 0.22695 0.2083 0.144192 0.15823 0.214212 0.171394 0.188113 0.173825 0.178753 0.199794 0.142583 0.187232 0.143863 0.361764 0.222843 0.186059 0.189758 0.224814 0.427179 0.208593 0.1568 0.277052 0.210007 0.17202 0.208786 0.326199 0.247512 0.204803 0.248991 0.338486 0.190687 0.182572 0.204072 0.262433 0.182311 0.166845 0.237176 0.174624 0.311511 0.218199 0.150481
Power_log_normal 0.194061 0.148582 0.201813 0.321874 0.194547 0.365181 0.153759 0.196309 0.160916 0.174986 0.164681 0.167747 0.266698 0.257267 0.16961 0.212385 0.143616 0.302767 0.180886 0.231085 0.159373 0.280174 0.221133 0.226794 0.216746 0.149719 0.215768 0.243774 0.162845 0.250382 0.217692 0.389052 0.246481 0.153778 0.294727 0.154745 0.356407 0.287725 0.294675 0.348688 0.331367 0.142264 0.304421 0.196408 0.226465 0.122179 0.283891 0.231626 0.272206 0.163052 0.236388 0.208418 0.334875 0.246287 0.245586 0.174214 0.225675 0.172893 0.252246 0.258415 0.289314 0.226972 0.37053 0.215846 0.166243 0.427872 0.318024 0.123228 0.245975 0.173654 0.164382 0.176308 0.377229 0.159805 0.194169 0.221017 0.24979 0.163991 0.253341 0.239065 0.387984 0.252064 0.179351 0.327613 0.254427 0.213276 0.142517 0.234736 0.280121 0.347075 0.229578 0.151847 0.206524 0.194102 0.247119 0.253406 0.408846 0.230282 0.272119 0.163039 0.254639 0.323491 0.212317 0.114684 0.113672 0.313035 0.170035 0.322682 0.202622 0.1903 0.186511 0.183414 0.243968 0.139322 0.398286 0.257037 0.49463 0.251017 0.195597 0.136178 0.208094 0.155213 0.178242 0.199471 0.230949 0.18851 0.267478 0.290869 0.186447 0.293622 0.254697 0.307197 0.365799 0.189313 0.244666 0.218899 0.147501 0.386799 0.211558 0.194207 0.179255 0.254611 0.187661 0.336258 0.259135 0.250079 0.172972 0.206352 0.241009 0.196877 0.208547 0.242479 0.341138 0.199548 0.301592 0.191011 0.154734 0.141069 0.360942 0.231556 0.236529 0.2153 0.313578 0.185713 0.267221 0.240351 0.316766 0.389578 0.259989 0.155227 0.188808 0.245041 0.266188 0.294321 0.303475 0.288924 0.233396 0.266408 0.180053 0.169155 0.146098 0.234987 0.28821 0.234936 0.31853 0.342468 0.266407 0.166464 0.31689 0.17406 0.301477 0.229922 0.173749 0.210134 0.131947 0.214039 0.218889 0.168488 0.331594 0.135431 0.190529 0.185166 0.178534 0.260564 0.268746 0.317893 0.124731 0.236842 0.258124 0.302792 0.175497 0.13432 0.205522 0.181897 0.253677 0.336307 0.255677 0.219683 0.1827 0.280443 0.127691 0.2107 0.221254 0.283669 0.232636 0.191457 0.243273 0.146438 0.208431 0.184162 0.169189 0.236556 0.241031 0.476151 0.172516 0.255847 0.376967 0.410678 0.298546 0.21785 0.1533 0.21671 0.293635 0.239331 0.127693 0.265372 0.265232 0.113341 0.277098 0.298095 0.109777 0.187735 0.267908 0.243636 0.136454 0.18172 0.283529 0.312918 0.442478 0.233946 0.163446 0.183973 0.306407 0.204805 0.242565 0.165359 0.295633 0.142289 0.296457 0.229512 0.18132 0.206546 0.281279 0.265742 0.216247 0.404215 0.176702 0.339142 0.299843 0.204674 0.192665 0.166249 0.230444 0.320496 0.174323 0.161618 0.269233 0.170467 0.255907 0.208998 0.279734 0.310381 0.185117 0.29642 0.162054 0.217888 0.229104 0.291779 0.193374 0.319339 0.249077 0.207422 0.243974 0.193129 0.119451 0.163299 0.218255 0.144519 0.162841 0.120914 0.267527 0.209084 0.200546 0.281479 0.273524 0.181079 0.179277 0.12943 0.256683 0.283308 0.202823 0.249496 0.235372 0.227666 0.224681 0.244024 0.168565 0.260349 0.260869 0.15364 0.224255 0.333648 0.136283 0.160133 0.193697 0.228614 0.19746 0.377259 0.174097 0.240412 0.205648 0.200587 0.181341 0.177235 0.26517 0.189412 0.201385 0.288244 0.285346 0.354194 0.21774 0.249675 0.209081 0.171833 0.173343 0.178397 0.259825 0.165829 0.301668 0.32817 0.216581 0.237363 0.228162 0.343358 0.142974 0.381257 0.187197 0.177907 0.196893 0.149284 0.212613 0.134969 0.103671 0.194524 0.304675 0.146428 0.208218 0.168932 0.314011 0.174447 0.223878 0.212264 0.258957 0.399415 0.258629 0.318199 0.26948 0.188819 0.176005 0.151841 0.218302 0.252121 0.286075 0.188087 0.154244 0.351006 0.2998 0.214005 0.252606 0.212905 0.242156 0.199143 0.24545 0.248748 0.31188 0.23646 0.208817 0.344747 0.176383 0.180726 0.195419 0.195209 0.226962 0.275821 0.152733 0.207295 0.15909 0.178098 0.213103 0.207935 0.444879 0.229296 0.42087 0.248339 0.232226 0.254625 0.260717 0.256295 0.148532 0.154827 0.226148 0.358467 0.212664 0.222281 0.140784 0.223698 0.33148 0.241284 0.371263 0.229047 0.22299 0.187172 0.342296 0.224435 0.194002 0.094232 0.387224 0.384703 0.192689 0.229279 0.149164 0.192191 0.220292 0.401259 0.341074 0.222611 0.155175 0.179461 0.143888 0.198041 0.23902 0.354506 0.247632 0.30736 0.217423 0.29269 0.281409 0.139321 0.251549 0.125312 0.120449 0.183917 0.181596 0.134148 0.167346 0.357408 0.308974 0.305638 0.135841 0.284674 0.315173 0.222293 0.158611 0.231781 0.166739 0.316087 0.298728 0.203217 0.222806 0.221703 0.219598 0.205839 0.28311 0.276392 0.224031 0.251402 0.212955 0.283797 0.335973 0.121833 0.221042 0.161115 0.176652 0.248377 0.363032 0.127371 0.19534 0.208788 0.205297 0.231496 0.273541 0.180768 0.179715 0.320434 0.223701 0.216772 0.283696 0.289676 0.215716 0.174993 0.293452 0.183965 0.149605 0.14762 0.195449 0.3464 0.225554 0.201983 0.304854 0.166936 0.200305 0.244495 0.219191 0.240322 0.151655 0.286419 0.19477 0.438791 0.136097 0.391173 0.325111 0.166247 0.241125 0.151147 0.221434 0.33414 0.213965 0.38309 0.184827 0.183367 0.259489 0.219377 0.179062 0.165603 0.201944 0.149974 0.172994 0.22312 0.191158 0.22469 0.369899 0.223964 0.174367 0.272793 0.262577 0.258497 0.139818 0.21554 0.338004 0.291351 0.289955 0.174621 0.345721 0.179771 0.268325 0.204847 0.243554 0.287656 0.121074 0.288864 0.199626 0.236031 0.167682 0.288541 0.178011 0.189142 0.186582 0.343191 0.216943 0.245474 0.145749 0.307118 0.199149 0.320622 0.155127 0.138837 0.306076 0.266826 0.335557 0.172083 0.24195 0.317989 0.165673 0.2795 0.206156 0.188237 0.235553 0.230993 0.250246 0.135356 0.360445 0.200127 0.165783 0.186452 0.236386 0.227868 0.191855 0.250999 0.138041 0.247817 0.20277 0.413117 0.306698 0.159118 0.315117 0.13596 0.167097 0.180877 0.15174 0.189902 0.253914 0.329004 0.236753 0.351084 0.17951 0.151445 0.167733 0.279727 0.209043 0.203613 0.302614 0.139546 0.240421 0.300046 0.209977 0.211168 0.41824 0.15996 0.345348 0.179112 0.190282 0.174394 0.274188 0.251394 0.306236 0.243712 0.37245 0.278405 0.295853 0.389492 0.257401 0.457908 0.486784 0.254272 0.289851 0.206337 0.212875 0.37814 0.202092 0.220523 0.171067 0.189986 0.197172 0.222225 0.107512 0.207132 0.387235 0.128599 0.188753 0.18943 0.23457 0.238976 0.199799 0.152454 0.258267 0.169454 0.254787 0.191757 0.217714 0.201316 0.240202 0.220429 0.209354 0.17008 0.336168 0.150137 0.153358 0.230238 0.248065 0.290377 0.319256 0.241267 0.193253 0.179759 0.153293 0.244706 0.270923 0.194648 0.225643 0.378008 0.29752 0.152494 0.206533 0.201215 0.317785 0.237021 0.197805 0.212571 0.203299 0.432294 0.275334 0.287912 0.232803 0.207604 0.173045 0.218368 0.094664 0.371856 0.246052 0.184684 0.306296 0.24365 0.236061 0.218105 0.209246 0.289694 0.396817 0.184299 0.204119 0.179664 0.257897 0.300498 0.162052 0.241301 0.237471 0.222512 0.213617 0.436912 0.202637 0.144443 0.22944 0.197366 0.405021 0.152537 0.221244 0.119254 0.220451 0.204424 0.27282 0.247614 0.220879 0.313957 0.34943 0.283831 0.261502 0.243909 0.324257 0.121596 0.142049 0.189451 0.166313 0.134444 0.200245 0.146421 0.269684 0.176084 0.244276 0.247852 0.226355 0.172283 0.213625 0.172587 0.255565 0.268416 0.218891 0.248411 0.194965 0.126602 0.182972 0.222212 0.240664 0.347446 0.283007 0.180287 0.236785 0.172611 0.209591 0.199458 0.380258 0.270802 0.149755 0.220599 0.142907 0.152471 0.358669 0.212669 0.180467 0.198395 0.282488 0.162456 0.210353 0.182372 0.178718 0.212784 0.289066 0.221939 0.242784 0.319047 0.225895 0.156687 0.129397 0.153132 0.252226 0.199745 0.278748 0.285294 0.211639 0.294984 0.195332 0.196427 0.399883 0.307412 0.305733 0.303474 0.295323 0.204017 0.460278 0.190232 0.207882 0.250544 0.210413 0.19788 0.240256 0.20999 0.278689 0.262021 0.304798 0.229122 0.197361 0.266111 0.17734 0.21058 0.355112 0.211827 0.173573 0.135111 0.469711 0.184062 0.184736 0.246418 0.322316 0.285088 0.21797 0.222491 0.261914 0.312019 0.210908 0.344013 0.339034 0.27262 0.251043 0.233167 0.324283 0.245784 0.176303 0.288642 0.268267 0.151869 0.177132 0.138374 0.21576 0.340795 0.219207 0.182531 0.207912 0.290647 0.17703 0.186217 0.285557 0.160251 0.18515 0.241688 0.334771 0.192685 0.148446 0.224299 0.218291 0.195431 0.335257 0.354764 0.227462 0.293689 0.256943 0.277458 0.311717 0.206714 0.363869 0.287532 0.196347 0.210217 0.178785 0.129498 0.193787 0.313983 0.204493 0.158387 0.150773 0.238189 0.29211 0.225929 0.276571 0.271386 0.138435 0.153451 0.238368 0.212275 0.213161 0.240147 0.150854 0.215347 0.164792 0.157848 0.203597 0.406662 0.306533 0.376956 0.223383 0.167609 0.153447 0.26878 0.244876 0.210462 0.184134 0.427906 0.193572 0.206604 0.233776 0.335494 0.249647 0.272171 0.33877 0.207261 0.216246 0.296838 0.221666 0.146121 0.211662 0.258816 0.195567 0.185993 0.322017 0.286617 0.188493 0.247432 0.255099 0.147697 0.14004 0.144009 0.137653 0.246093 0.274681 0.139816 0.247558 0.210254 0.158516 0.163744 0.179708 0.20892 0.229893 0.322932 0.213428 0.225516 0.253109 0.270686 0.431068 0.292925 0.192449 0.197903 0.199049 0.148589 0.177669 0.247914 0.228512 0.204925 0.114586 0.172823 0.365427 0.260718 0.186177 0.169986 0.27362 0.1163 0.285731 0.282853 0.216886 0.192662 0.161858 0.231687 0.274991 0.251015 0.289454 0.184755 0.202197 0.259383 0.183213 0.173369 0.208794 0.230214 0.20044 0.385044 0.248017 0.233681 0.204665 0.303722 0.219652 0.351823 0.343129 0.434725 0.332499 0.237335 0.174562 0.192296 0.240541 0.177896 0.268992 0.185374 0.273958 0.243739 0.273499 0.232401 0.189943 0.335249 0.267158 0.172758 0.143643 0.198161 0.225891 0.213082 0.253335 0.181435 0.397256 0.240487 0.129015 0.219642 0.261601 0.378479 0.226634 0.223935 0.296401 0.168252 0.199207 0.195618 0.226644 0.349336 0.121874 0.292039 0.355483 0.254864 0.328881 0.273697 0.251269 0.13558 0.119721 0.38278 0.217763 0.260369 0.123932 0.226999 0.234151 0.204861 0.143283 0.271852 0.232756 0.381551 0.20486 0.256312 0.277443 0.172501 0.177384 0.237568 0.340522 0.143556 0.235033 0.193763 0.148837 0.186553 0.184496 0.214569 0.348265 0.231304 0.220395 0.200586 0.22559 0.252771 0.216465 0.34878 0.217686 0.149153 0.351075 0.210321 0.240034 0.242964 0.293132 0.23609 0.219558 0.219206 0.271639 0.275317 0.120591 0.395665 0.227739 0.145404 0.226583 0.217102 0.287572 0.133564 0.203332 0.224693 0.359784 0.285232 0.298864 0.122965 0.231819 0.202566 0.333872 0.290825 0.231789 0.243203 0.235783 0.194603 0.291761 0.171309 0.163852 0.14163 0.210037 0.186444 0.299915 0.296877 0.102197 0.251731 0.246155 0.214287 0.194937 0.147319 0.330204 0.180855 0.211924 0.171268 0.175132 0.371181 0.441652 0.294142 0.26153 0.296653 0.163031 0.160082 0.252974 0.26989 0.150284 0.256991 0.172175 0.239996 0.286153 0.289138 0.113261 0.227588 0.255071 0.158422 0.426609 0.225683 0.250048 0.192269 0.196858 0.178682 0.146993 0.246023 0.192909 0.20653 0.264739 0.268332 0.202331 0.205658 0.315087 0.179014 0.360189 0.190556 0.240902 0.192159 0.200174 0.260915 0.402233 0.266174 0.170486 0.131021 0.491612 0.19426 0.203652 0.273985 0.222928 0.206497 0.153949 0.347695 0.14494 0.231758 0.204738 0.206404 0.212846 0.19149 0.223575 0.367849 0.186796 0.233438 0.154733 0.172021 0.298369 0.176255 0.38203 0.257875 0.243675 0.233556 0.273154 0.134228 0.225742 0.344473 0.348228 0.293069 0.252813 0.159347 0.141112 0.218744 0.15389 0.230512 0.235556 0.20651 0.279552 0.194359 0.327487 0.146395 0.222454 0.377634 0.159255 0.230964 0.339922 0.22846 0.370852 0.203645 0.226299 0.289986 0.28323 0.183821 0.298404 0.239393 0.185775 0.175776 0.1908 0.1985 0.279518 0.256793 0.283909 0.243201 0.210236 0.144486 0.311975 0.250646 0.277306 0.132691 0.281522 0.277632 0.31994 0.345083 0.194705 0.204138 0.353861 0.247651 0.259791 0.160832 0.351834 0.309917 0.173762 0.207453 0.229404 0.252083 0.232004 0.176411 0.195747 0.258105 0.239105 0.21068 0.161789 0.171623 0.1503 0.287172 0.391738 0.229061 0.271124 0.243765 0.189655 0.15432 0.12262 0.206827 0.217755 0.224828 0.263183 0.361243 0.315241 0.16873 0.193534 0.174226 0.191131 0.205075 0.107273 0.222372 0.199819 0.160494 0.132693 0.213414 0.243194 0.325651 0.286714 0.298234 0.192133 0.286678 0.162284 0.219 0.316736 0.179621 0.130822 0.223232 0.446585 0.19722 0.236223 0.235199 0.213924 0.230838 0.167428 0.177842 0.221539 0.295337 0.477571 0.233941 0.19307 0.165683 0.140532 0.264413 0.27906 0.20198 0.225923 0.357967 0.319858 0.270884 0.247194 0.394048 0.355192 0.214967 0.231218 0.227028 0.126372 0.349967 0.236028 0.263276 0.243881 0.243911 0.22563 0.204051 0.180543 0.245436 0.240863 0.235533 0.16512 0.209478 0.296873 0.119556 0.229547 0.15654 0.232206 0.135109 0.264349 0.16306 0.254592 0.24534 0.221918 0.189239 0.294421 0.158444 0.21436 0.249495 0.241417 0.277898 0.3741 0.150497 0.201632 0.220046 0.306793 0.306292 0.235469 0.104635 0.214822 0.274256 0.134549 0.202517 0.182144 0.220649 0.311759 0.29001 0.291962 0.273107 0.296156 0.301682 0.207327 0.181453 0.221205 0.29213 0.200574 0.230893 0.185124 0.255201 0.151398 0.213143 0.152114 0.202053 0.213756 0.206552 0.195396 0.318074 0.221788 0.148616 0.208019 0.346822 0.231028 0.136888 0.139421 0.225212 0.253778 0.185159 0.12696 0.132779 0.177202 0.307867 0.182935 0.112674 0.428511 0.19281 0.181515 0.246976 0.191501 0.232815 0.209801 0.218331 0.223822 0.223168 0.273846 0.193908 0.324304 0.303156 0.194894 0.385384 0.271626 0.339872 0.229738 0.269175 0.242641 0.371982 0.291869 0.105257 0.270349 0.188093 0.215394 0.361019 0.224093 0.166007 0.193895 0.190443 0.139917 0.249397 0.186204 0.306852 0.314066 0.160147 0.164149 0.365337 0.128363 0.218266 0.244884 0.219878 0.148056 0.19224 0.252346 0.1602 0.215217 0.18785 0.225555 0.120312 0.196064 0.233102 0.310852 0.309821 0.206371 0.188093 0.184277 0.21848 0.186762 0.154557 0.321238 0.243165 0.321537 0.234036 0.20849 0.199602 0.233779 0.198353 0.245883 0.220302 0.236589 0.195224 0.18167 0.242238 0.178346 0.193374 0.242643 0.233762 0.228304 0.232946 0.304481 0.383801 0.214295 0.28283 0.239989 0.26887 0.239162 0.373408 0.244786 0.134468 0.160699 0.264315 0.204325 0.263201 0.304998 0.229305 0.381726 0.154393 0.234432 0.156929 0.194488 0.152737 0.144832 0.288308 0.216546 0.336281 0.297364 0.144533 0.230315 0.201423 0.286197 0.345562 0.156439 0.219873 0.377228 0.201571 0.286446 0.143382 0.143411 0.2759 0.19962 0.220037 0.169287 0.171971 0.267825 0.155207 0.302667 0.248571 0.247156 0.24922 0.293325 0.314946 0.136046 0.215643 0.282065 0.093702 0.222594 0.316267 0.225973 0.198232 0.152877 0.207061 0.359306 0.097488 0.255444 0.202276 0.211613 0.391969 0.211963 0.25391 0.130615 0.271493 0.288604 0.216904 0.177573 0.23506 0.204466 0.270943 0.409588 0.270099 0.198869 0.188646 0.288936 0.215642 0.376279 0.163968 0.357935 0.152255 0.192075 0.362055 0.178644 0.395847 0.159139 0.163639 0.343102 0.300527 0.329444 0.286495 0.331446 0.284237 0.331254 0.242597 0.25247 0.159994 0.271857 0.212715 0.238276 0.266222 0.250815 0.247156 0.354081 0.175977 0.167704 0.216565 0.262246 0.226988 0.228278 0.35712 0.241397 0.273316 0.257553 0.21623 0.25562 0.292086 0.205573 0.25534 0.236151 0.158256 0.14534 0.173072 0.33251 0.40125 0.15892 0.318955 0.206496 0.246587 0.256146 0.194778 0.263937 0.194919 0.27859 0.232159 0.156689 0.28095 0.218564 0.175433 0.386384 0.177994 0.22147 0.243872 0.15722 0.173782 0.141022 0.225054 0.157175 0.329657 0.219753 0.223871 0.173225 0.180999 0.221499 0.191494 0.21993 0.238977 0.144842 0.293158 0.341759 0.234639 0.300667 0.213906 0.287653 0.220532 0.227553 0.233611 0.160416 0.231738 0.256693 0.137315 0.25893 0.217616 0.216715 0.231691 0.112864 0.254491 0.295977 0.30411 0.354332 0.234076 0.211897 0.308994 0.168837 0.154965 0.215858 0.243267 0.237387 0.18931 0.195057 0.233754 0.244736 0.273369 0.184415 0.275156 0.224994 0.229258 0.27516 0.145328 0.232914 0.278162 0.217205 0.210054 0.153383 0.447651 0.135503 0.239549 0.229155 0.196228 0.196302 0.209432 0.210394 0.25678 0.190261 0.18511 0.249239 0.243268 0.200108 0.205311 0.113749 0.256062 0.23525 0.158466 0.332637 0.240163 0.154737 0.21629 0.138477 0.281308 0.168011 0.273533 0.203106 0.135718 0.212705 0.261747 0.16825 0.212108 0.282769 0.366388 0.230004 0.191254 0.138376 0.140797 0.384984 0.217274 0.242534 0.205447 0.255115 0.214848 0.216511 0.156548 0.204879 0.471354 0.208121 0.177793 0.282083 0.227344 0.133897 0.225162 0.131651 0.34488 0.264686 0.309934 0.325911 0.161221 0.223399 0.199154 0.253366 0.202555 0.323919 0.210755 0.240141 0.181304 0.321059 0.192522 0.313128 0.311773 0.205055 0.387368 0.300792 0.14181 0.205519 0.277379 0.269861 0.147546 0.216745 0.179004 0.28244 0.181064 0.145925 0.20613 0.169911 0.317366 0.128621 0.318415 0.165586 0.259014 0.203849 0.14234 0.259781 0.250842 0.18802 0.259817 0.347714 0.179801 0.192918 0.215523 0.320697 0.274143 0.249453 0.18339 0.272707 0.194577 0.192664 0.294783 0.20227 0.290382 0.410867 0.221647 0.27877 0.178197 0.354278 0.270425 0.27265 0.24796 0.204153 0.210963 0.150694 0.313914 0.241086 0.289133 0.244384 0.172322 0.315956 0.287149 0.236395 0.279735 0.135857 0.24891 0.244435 0.21598 0.228835 0.209713 0.145787 0.181608 0.396209 0.142706 0.196809 0.266778 0.184167 0.213672 0.179537 0.310482 0.329304 0.234443 0.302117 0.209299 0.336169 0.201633 0.276137 0.192777 0.271607 0.510602 0.351558 0.529858 0.403635 0.194731 0.219112 0.151737 0.178688 0.241068 0.295644 0.168408 0.153495 0.284988 0.163894 0.274272 0.176667 0.207815 0.185459 0.266837 0.308258 0.220756 0.190298 0.207036 0.262969 0.191601 0.252232 0.266365 0.215202 0.345336 0.225689 0.173591 0.222643 0.120297 0.311098 0.202045 0.256079 0.277308 0.165748 0.197727 0.195249 0.142771 0.179942 0.139369 0.25898 0.270639 0.265565 0.365606 0.342957 0.191297 0.193893 0.143866 0.176 0.288321 0.121486 0.171203 0.126302 0.237144 0.16913 0.139258 0.134354 0.25255 0.26165 0.370334 0.259942 0.200867 0.210596 0.215327 0.299293 0.299521 0.215289 0.202498 0.160701 0.195077 0.163613 0.231899 0.182087 0.280433 0.1795 0.174524 0.269347 0.305894 0.259943 0.118407 0.260749 0.291401 0.211475 0.157538 0.214387 0.137134 0.160786 0.252195 0.271529 0.195047 0.217713 0.321805 0.323448 0.204191 0.311793 0.271088 0.326435 0.09261 0.135436 0.20023 0.381852 0.146032 0.202671 0.132588 0.294089 0.322459 0.195941 0.240303 0.19945 0.265231 0.180621 0.184898 0.163871 0.160839 0.326603 0.205703 0.139391 0.177792 0.162593 0.24633 0.161794 0.240624 0.32977 0.223703 0.13961 0.256988 0.247042 0.267061 0.381233 0.126673 0.28286 0.225016 0.214342 0.423079 0.196736 0.252718 0.229282 0.270865 0.321821 0.261827 0.174843 0.192061 0.16943 0.220584 0.213511 0.179664 0.304993 0.205356 0.14827 0.219658 0.163032 0.204856 0.305957 0.172224 0.223414 0.286351 0.19154 0.451487 0.169739 0.228316 0.177841 0.232318 0.238073 0.216146 0.148873 0.116383 0.218145 0.136062 0.222003 0.198495 0.180494 0.318615 0.161579 0.258057 0.184301 0.293518 0.217305 0.288565 0.26234 0.264267 0.318377 0.250695 0.223703 0.240992 0.283437 0.171171 0.160723 0.225988 0.287474 0.185643 0.324862 0.216508 0.248147 0.149215 0.216894 0.194744 0.159037 0.177638 0.27478 0.173148 0.172313 0.343343 0.221408 0.202693 0.188675 0.204078 0.114762 0.154793 0.217262 0.13482 0.315008 0.216444 0.208956 0.175925 0.25036 0.125034 0.101022 0.158965 0.184901 0.153211 0.349828 0.422778 0.155379 0.166492 0.15623 0.211264 0.280833 0.251589 0.331436 0.334979 0.261377 0.179989 0.191906 0.442604 0.25569 0.214895 0.249535 0.238888 0.207736 0.167433 0.264896 0.161872 0.2985 0.258754 0.158124 0.282663 0.254366 0.253085 0.189302 0.182726 0.143198 0.170322 0.326073 0.20958 0.241542 0.264023 0.236017 0.175485 0.077725 0.234034 0.219043 0.218524 0.20699 0.160897 0.221045 0.252561 0.410563 0.346624 0.350666 0.228455 0.208533 0.178872 0.380967 0.213139 0.166988 0.258037 0.207783 0.293206 0.250224 0.19723 0.207033 0.218464 0.218887 0.302256 0.270516 0.240637 0.306898 0.113453 0.151689 0.22 0.316374 0.157245 0.185564 0.191807 0.156788 0.220491 0.35543 0.23287 0.213191 0.239682 0.161229 0.272696 0.163628 0.159439 0.285394 0.23301 0.252262 0.199658 0.241851 0.166511 0.182623 0.268425 0.287277 0.159766 0.308475 0.30217 0.139377 0.29153 0.307791 0.223442 0.149223 0.262422 0.21043 0.136258 0.205152 0.237814 0.14845 0.229323 0.263278 0.233979 0.276892 0.188104 0.290023 0.301987 0.200192 0.262642 0.239411 0.158046 0.277969 0.200664 0.223587 0.155968 0.146521 0.149369 0.296482 0.272554 0.166347 0.185363 0.24345 0.381753 0.139085 0.194941 0.171963 0.21701 0.146038 0.333155 0.150925 0.164363 0.242956 0.222906 0.2886 0.173867 0.264894 0.242969 0.314148 0.229997 0.20998 0.306195 0.165697 0.210911 0.211811 0.161223 0.233282 0.286703 0.167627 0.163774 0.257356 0.220667 0.357566 0.121686 0.270547 0.237373 0.208535 0.202622 0.249264 0.202035 0.160755 0.244032 0.284066 0.245629 0.124001 0.269421 0.217407 0.249095 0.244267 0.225912 0.30843 0.332944 0.325237 0.179807 0.168263 0.477001 0.206967 0.265182 0.26859 0.230285 0.23031 0.288193 0.204461 0.295887 0.253149 0.173651 0.226746 0.292268 0.205831 0.283852 0.365325 0.192251 0.261257 0.310868 0.329416 0.254367 0.110091 0.376069 0.171039 0.157253 0.2912 0.272717 0.228546 0.377694 0.290323 0.194184 0.241157 0.184599 0.43332 0.170003 0.353358 0.251299 0.301511 0.19284 0.212491 0.218915 0.227427 0.17293 0.265761 0.201496 0.208804 0.274872 0.278303 0.287403 0.287272 0.301409 0.214605 0.167399 0.228814 0.200316 0.21602 0.327392 0.381918 0.214686 0.240846 0.188145 0.165911 0.343432 0.194472 0.197801 0.204578 0.245544 0.277742 0.250834 0.276628 0.211206 0.209226 0.193583 0.268739 0.350034 0.23646 0.23223 0.228537 0.269869 0.235682 0.314826 0.47205 0.213592 0.346295 0.272931 0.296525 0.193706 0.145049 0.212541 0.182567 0.245303 0.154485 0.237751 0.343547 0.226108 0.32084 0.285755 0.219305 0.137052 0.238058 0.161261 0.276087 0.336552 0.299293 0.31875 0.273613 0.227671 0.339882 0.237961 0.335406 0.267079 0.240622 0.26408 0.255927 0.236463 0.179597 0.255997 0.183695 0.171263 0.468005 0.336614 0.190563 0.229108 0.210159 0.294334 0.167712 0.317222 0.247596 0.207777 0.217274 0.318089 0.209818 0.326526 0.188183 0.218085 0.186877 0.327399 0.179168 0.187867 0.300867 0.284196 0.207775 0.25809 0.179973 0.283144 0.180219 0.286649 0.250239 0.203955 0.262094 0.272714 0.266335 0.313561 0.156763 0.254234 0.240558 0.227728 0.293304 0.104142 0.261149 0.131553 0.273315 0.152826 0.500079 0.194087 0.286315 0.241018 0.208749 0.355389 0.157211 0.263038 0.242682 0.347788 0.157709 0.223906 0.203478 0.334358 0.22007 0.249529 0.292708 0.175103 0.228886 0.183337 0.182111 0.256015 0.181361 0.233203 0.323783 0.175595 0.217176 0.175191 0.140715 0.185705 0.283362 0.192798 0.187796 0.160762 0.241525 0.253211 0.231707 0.112663 0.180384 0.338097 0.115938 0.261825 0.336918 0.321992 0.259704 0.237951 0.21383 0.282723 0.139436 0.185705 0.126983 0.205261 0.265673 0.106159 0.191276 0.102776 0.26048 0.184994 0.261363 0.188531 0.211368 0.208414 0.221936 0.267537 0.266996 0.236197 0.229358 0.163856 0.232782 0.189143 0.167065 0.153031 0.178027 0.26035 0.176492 0.270558 0.169377 0.180597 0.434372 0.246504 0.306967 0.207202 0.242274 0.368103 0.3465 0.122837 0.240727 0.226262 0.243162 0.374817 0.338314 0.200386 0.267439 0.219986 0.21613 0.319612 0.143463 0.212134 0.167699 0.113017 0.14837 0.221026 0.244916 0.144241 0.242033 0.17249 0.216503 0.173317 0.278729 0.185302 0.22902 0.245941 0.240241 0.214791 0.202614 0.273521 0.188356 0.289746 0.208315 0.131413 0.196071 0.234198 0.2864 0.249649 0.247687 0.177602 0.182996 0.234998 0.135912 0.330439 0.287304 0.162969 0.32453 0.250671 0.225903 0.407558 0.255409 0.265436 0.163502 0.245923 0.215427 0.219118 0.277349 0.198255 0.154057 0.267771 0.248849 0.185057 0.183039 0.280945 0.196214 0.337619 0.148232 0.338693 0.191445 0.211371 0.166027 0.197047 0.164592 0.221476 0.262078 0.305475 0.330967 0.174168 0.206472 0.223407 0.510937 0.209283 0.219232 0.332522 0.156812 0.256039 0.219633 0.236593 0.419081 0.313795 0.218688 0.220311 0.287862 0.150032 0.166559 0.371932 0.11896 0.272571 0.396748 0.349642 0.226071 0.246496 0.207269 0.348144 0.298504 0.249319 0.268779 0.22305 0.373783 0.222939 0.214806 0.161886 0.272054 0.227662 0.350835 0.246754 0.270249 0.31609 0.181429 0.178031 0.119023 0.112389 0.161554 0.279634 0.232469 0.26415 0.219015 0.188728 0.254598 0.28635 0.2209 0.194587 0.255999 0.257405 0.243193 0.225875 0.183232 0.353356 0.167105 0.248687 0.179348 0.345044 0.407169 0.29169 0.15407 0.164789 0.206747 0.27634 0.366301 0.273847 0.238502 0.15687 0.217268 0.190028 0.17601 0.310734 0.181575 0.247234 0.283654 0.246248 0.193187 0.294835 0.157377 0.262513 0.257476 0.292306 0.242479 0.22233 0.2188 0.17427 0.294639 0.209147 0.155017 0.175411 0.310189 0.077481 0.306607 0.201895 0.347764 0.185579 0.283162 0.212318 0.218406 0.233605 0.246993 0.205042 0.176419 0.15321 0.224684 0.341125 0.17438 0.219028 0.346552 0.332843 0.258278 0.272541 0.226574 0.380491 0.40243 0.140548 0.202154 0.209095 0.249199 0.174738 0.272024 0.319537 0.196417 0.248049 0.313512 0.257195 0.46729 0.338115 0.175556 0.176038 0.206059 0.275426 0.36124 0.275605 0.218301 0.16326 0.222391 0.302062 0.291115 0.212931 0.373817 0.16584 0.143777 0.252821 0.228902 0.107673 0.274395 0.314434 0.256583 0.316957 0.203709 0.174135 0.188856 0.347062 0.30139 0.298742 0.25369 0.20336 0.11503 0.317813 0.18512 0.200259 0.175072 0.374877 0.230266 0.270274 0.164056 0.2178 0.256214 0.313657 0.247759 0.227505 0.252559 0.162669 0.219314 0.229019 0.302695 0.392249 0.246703 0.381231 0.191815 0.354034 0.293313 0.303591 0.127383 0.227088 0.154558 0.190876 0.18902 0.257995 0.196225 0.28145 0.30882 0.286389 0.160966 0.394108 0.281084 0.145307 0.156401 0.328022 0.174002 0.290483 0.191595 0.196248 0.205015 0.430856 0.179206 0.172012 0.160009 0.261543 0.199234 0.256974 0.16167 0.253471 0.287036 0.184106 0.219539 0.236844 0.226182 0.248193 0.343742 0.254926 0.216522 0.382063 0.207724 0.143735 0.266346 0.159609 0.301862 0.155958 0.173965 0.236569 0.173915 0.216079 0.236128 0.223033 0.220763 0.173709 0.218459 0.302819 0.256608 0.205468 0.184826 0.205508 0.218495 0.131355 0.258377 0.420274 0.161122 0.174592 0.227674 0.272445 0.186622 0.311186 0.316428 0.132176 0.188809 0.157628 0.308941 0.230032 0.279075 0.165372 0.226518 0.287372 0.344026 0.143118 0.232096 0.299558 0.229322 0.258173 0.226404 0.155059 0.166831 0.156294 0.421314 0.197486 0.234615 0.24113 0.177909 0.333836 0.305298 0.140824 0.250163 0.191272 0.153167 0.311364 0.227152 0.359199 0.234988 0.147065 0.301907 0.319724 0.282261 0.125455 0.206988 0.229675 0.112837 0.249829 0.221144 0.197455 0.190755 0.215052 0.185452 0.203566 0.228916 0.335276 0.189421 0.261045 0.28176 0.188339 0.148853 0.148779 0.378427 0.359535 0.35467 0.164541 0.140386 0.129883 0.139241 0.131282 0.256983 0.139015 0.145203 0.331869 0.128234 0.161869 0.298413 0.178775 0.238921 0.283215 0.222622 0.303434 0.198 0.183824 0.223187 0.22244 0.237082 0.396575 0.187746 0.222219 0.167496 0.210039 0.253149 0.138648 0.380111 0.17672 0.139641 0.21705 0.205015 0.208504 0.223827 0.162273 0.269704 0.233514 0.135679 0.184945 0.147509 0.285942 0.217947 0.19081 0.189366 0.414438 0.268583 0.314933 0.227502 0.115702 0.202349 0.169374 0.317381 0.185154 0.161551 0.278324 0.2009 0.19412 0.183522 0.110321 0.247749 0.244478 0.130103 0.300039 0.297452 0.340134 0.129735 0.289749 0.172981 0.284649 0.296963 0.358956 0.19812 0.200089 0.230552 0.221859 0.255092 0.213253 0.18769 0.129136 0.23808 0.195884 0.146961 0.302751 0.162965 0.267849 0.203668 0.204944 0.151787 0.26501 0.165712 0.196501 0.286546 0.213733 0.128836 0.379673 0.204215 0.222794 0.236779 0.17239 0.213536 0.186037 0.212225 0.335488 0.372603 0.166903 0.268817 0.315 0.234786 0.268026 0.370845 0.309458 0.210102 0.276037 0.184636 0.218383 0.192733 0.104215 0.206572 0.178733 0.321916 0.283927 0.232435 0.22193 0.286198 0.253842 0.350909 0.255565 0.164753 0.250753 0.210919 0.38006 0.180862 0.362234 0.296737 0.370948 0.317999 0.233747 0.194886 0.220686 0.237901 0.165467 0.233558 0.123458 0.192301 0.158937 0.293387 0.342531 0.224277 0.179922 0.186348 0.116159 0.281468 0.20089 0.382892 0.180809 0.247525 0.321955 0.234291 0.247951 0.301264 0.140746 0.145733 0.203537 0.208633 0.27916 0.238168 0.139353 0.171812 0.395858 0.265151 0.191965 0.320858 0.181085 0.203125 0.240438 0.180506 0.257193 0.238504 0.415835 0.221349 0.204088 0.219813 0.182299 0.140134 0.283565 0.279686 0.242118 0.383555 0.338647 0.209832 0.277492 0.22836 0.244432 0.178864 0.369685 0.37293 0.215594 0.221793 0.212373 0.19016 0.377289 0.211781 0.310647 0.212272 0.286056 0.300506 0.28637 0.365521 0.136985 0.269119 0.16077 0.16425 0.383492 0.252422 0.386532 0.187894 0.142533 0.219395 0.257234 0.389674 0.312311 0.2303 0.247817 0.168755 0.331285 0.383583 0.141035 0.245186 0.258698 0.159999 0.191659 0.230687 0.147241 0.183451 0.19448 0.281007 0.185348 0.210372 0.280525 0.163842 0.206658 0.230001 0.229053 0.247956 0.145302 0.143519 0.223638 0.177484 0.245141 0.390638 0.204685 0.382854 0.299277 0.248052 0.215327 0.260416 0.182734 0.232836 0.30621 0.188957 0.225249 0.199446 0.223405 0.17838 0.304295 0.30617 0.270955 0.146271 0.134566 0.192984 0.200533 0.303742 0.290688 0.124435 0.323131 0.211662 0.299213 0.3031 0.168241 0.239487 0.219431 0.173901 0.356088 0.206165 0.354704 0.177568 0.131354 0.183366 0.312126 0.206242 0.144523 0.193653 0.244235 0.190732 0.139655 0.240926 0.255444 0.186018 0.311517 0.350577 0.196929 0.229284 0.296542 0.178062 0.169515 0.089211 0.233091 0.231703 0.352186 0.282426 0.246436 0.207772 0.163602 0.220659 0.289558 0.20967 0.146196 0.373379 0.278402 0.26024 0.167747 0.522037 0.218371 0.242344 0.275501 0.193092 0.251474 0.203657 0.183825 0.232897 0.115671 0.310716 0.193837 0.156416 0.396677 0.234646 0.139194 0.245015 0.152879 0.197549 0.351939 0.184438 0.168167 0.287234 0.148492 0.190878 0.2932 0.28893 0.276028 0.198864 0.396812 0.246405 0.17593 0.242011 0.33952 0.265293 0.197826 0.109526 0.360629 0.237764 0.259329 0.276955 0.298405 0.259309 0.221894 0.227154 0.274026 0.277689 0.262528 0.173881 0.201548 0.275103 0.203081 0.304023 0.203117 0.22523 0.162029 0.205759 0.377753 0.265191 0.211533 0.309625 0.254447 0.224002 0.423137 0.221699 0.168372 0.235001 0.161215 0.242943 0.17106 0.303173 0.221941 0.256642 0.363277 0.159095 0.234276 0.237499 0.187649 0.217853 0.259975 0.307135 0.200407 0.250001 0.2565 0.270969 0.254424 0.203019 0.344456 0.286418 0.303533 0.189041 0.277399 0.234148 0.46122 0.230747 0.180069 0.187823 0.259449 0.26759 0.217797 0.242692 0.360468 0.174045 0.248285 0.381902 0.211626 0.229642 0.23735 0.29038 0.347104 0.241456 0.249635 0.303231 0.219373 0.217466 0.318294 0.177764 0.201459 0.368423 0.138585 0.21591 0.215159 0.170129 0.212076 0.202989 0.158317 0.197863 0.234314 0.22162 0.359336 0.244601 0.21503 0.25798 0.312382 0.237407 0.307455 0.369416 0.118382 0.241206 0.184906 0.22402 0.206689 0.377708 0.303264 0.151549 0.270506 0.263921 0.187899 0.141951 0.296647 0.161787 0.194073 0.285144 0.202053 0.415034 0.270147 0.238353 0.264648 0.194076 0.303331 0.198172 0.339754 0.153181 0.137873 0.154528 0.19262 0.275888 0.160333 0.23314 0.382805 0.095097 0.148814 0.432323 0.321538 0.455797 0.258835 0.254401 0.296966 0.206873 0.277081 0.259867 0.167735 0.162674 0.232072 0.266368 0.149914 0.285176 0.241706 0.204401 0.220061 0.282122 0.187043 0.293462 0.205916 0.14926 0.211655 0.230363 0.181984 0.363583 0.222445 0.189662 0.185415 0.361432 0.193075 0.305476 0.181667 0.18409 0.153105 0.207703 0.316304 0.137832 0.250773 0.299988 0.226335 0.233381 0.229356 0.207482 0.18979 0.28721 0.299423 0.359301 0.292107 0.239937 0.248666 0.201457 0.154695 0.248545 0.181983 0.204418 0.376729 0.22758 0.350612 0.220824 0.191977 0.254452 0.274204 0.388805 0.190751 0.218522 0.241041 0.210832 0.279361 0.144204 0.240245 0.276519 0.241009 0.206356 0.419418 0.276373 0.312056 0.125713 0.214075 0.249213 0.511791 0.284976 0.31438 0.253084 0.134674 0.130478 0.201481 0.293634 0.293905 0.208037 0.179468 0.213828 0.142489 0.201783 0.175592 0.297916 0.17629 0.284417 0.200805 0.338953 0.262104 0.262205 0.193156 0.20044 0.366407 0.21526 0.276688 0.243833 0.181641 0.198728 0.174835 0.202163 0.389335 0.14275 0.194229 0.410842 0.299029 0.242555 0.323605 0.360056 0.248988 0.167474 0.251282 0.151046 0.253284 0.20882 0.272238 0.216001 0.206276 0.211862 0.394592 0.191383 0.149903 0.170392 0.200261 0.316162 0.176411 0.320199 0.208883 0.183213 0.15889 0.230173 0.145405 0.257994 0.340349 0.270537 0.335692 0.318128 0.140369 0.153365 0.256759 0.278348 0.258834 0.260685 0.165543 0.117151 0.179838 0.156068 0.189865 0.235793 0.26006 0.288579 0.246109 0.157593 0.185655 0.293178 0.243138 0.152433 0.197933 0.317423 0.231902 0.25731 0.203017 0.167131 0.182616 0.140983 0.293063 0.149728 0.277108 0.269752 0.167223 0.221471 0.27979 0.152906 0.100185 0.279286 0.186357 0.276289 0.289869 0.188975 0.375334 0.255077 0.34452 0.210642 0.187373 0.295805 0.218979 0.159247 0.20643 0.246633 0.265055 0.287657 0.203718 0.135194 0.366004 0.244647 0.246983 0.13048 0.158285 0.189108 0.260343 0.219422 0.26947 0.203826 0.375458 0.162492 0.195931 0.298163 0.153641 0.293828 0.393708 0.196017 0.1811 0.256452 0.261208 0.17641 0.121332 0.20062 0.283962 0.219207 0.1768 0.179699 0.206292 0.238319 0.215439 0.154287 0.202054 0.338635 0.224132 0.112691 0.241999 0.281342 0.227423 0.277939 0.254813 0.274259 0.300488 0.20573 0.243767 0.399703 0.181797 0.287642 0.230461 0.344083 0.129258 0.245952 0.121271 0.120233 0.303079 0.173084 0.190023 0.221584 0.233438 0.23096 0.35511 0.312641 0.171758 0.185189 0.137257 0.247704 0.285135 0.149649 0.338143 0.385114 0.19413 0.176461 0.21079 0.343752 0.1991 0.164614 0.241783 0.272222 0.172661 0.099639 0.295855 0.193843 0.26131 0.264672 0.196798 0.283211 0.245635 0.279456 0.305512 0.41271 0.169129 0.159395 0.179598 0.185979 0.251349 0.241974 0.171747 0.27704 0.350673 0.150695 0.195634 0.256243 0.273734 0.13848 0.286936 0.230691 0.164363 0.140324 0.339081 0.253571 0.273743 0.180865 0.247371 0.234605 0.38295 0.336631 0.189949 0.296624 0.23494 0.32117 0.427727 0.183123 0.194008 0.151963 0.170849 0.215491 0.173088 0.16983 0.204236 0.280078 0.18166 0.128141 0.328809 0.304976 0.233665 0.255213 0.186216 0.26727 0.262788 0.189126 0.187786 0.194505 0.253721 0.217262 0.2751 0.244605 0.140031 0.121991 0.170822 0.226806 0.301672 0.208795 0.187782 0.21569 0.21886 0.242696 0.201968 0.166768 0.225408 0.114409 0.25442 0.197569 0.231157 0.26983 0.166482 0.248005 0.234278 0.249412 0.346155 0.186211 0.217241 0.256029 0.163305 0.277677 0.144031 0.160438 0.29701 0.24509 0.21774 0.277497 0.232041 0.202757 0.18642 0.117713 0.225215 0.318926 0.25876 0.171154 0.184458 0.246412 0.207582 0.262625 0.351177 0.197339 0.309556 0.129462 0.179184 0.142607 0.337468 0.179382 0.230749 0.17465 0.169292 0.220349 0.12648 0.237917 0.17002 0.25239 0.382399 0.150823 0.271465 0.244027 0.340895 0.196874 0.30603 0.272603 0.158962 0.228488 0.246298 0.304493 0.244036 0.283962 0.325234 0.237224 0.27872 0.392193 0.224783 0.225905 0.199275 0.203404 0.209255 0.238169 0.280889 0.293253 0.217594 0.194291 0.305883 0.147349 0.261379 0.301656 0.251466 0.243696 0.412489 0.215084 0.255032 0.25002 0.146996 0.144847 0.222852 0.326035 0.173562 0.29704 0.19197 0.130822 0.269892 0.307441 0.145492 0.193247 0.176608 0.231622 0.349578 0.154216 0.227925 0.206834 0.144056 0.189955 0.181838 0.342237 0.242411 0.213838 0.198759 0.21417 0.339587 0.367009 0.144179 0.176468 0.180067 0.267391 0.230319 0.188629 0.213215 0.31492 0.240497 0.228096 0.257304 0.308595 0.372607 0.275493 0.217375 0.325781 0.289226 0.345784 0.325194 0.272354 0.229682 0.151768 0.160209 0.33488 0.243253 0.182132 0.154058 0.222295 0.303877 0.252749 0.185087 0.202218 0.220299 0.214017 0.42892 0.172736 0.219608 0.29075 0.23146 0.180363 0.265571 0.139221 0.208582 0.117523 0.237173 0.32983 0.259106 0.136492 0.286942 0.303255 0.27607 0.322749 0.181067 0.199112 0.238401 0.297421 0.221387 0.251063 0.136302 0.44994 0.176579 0.546447 0.29792 0.287864 0.118769 0.360331 0.181972 0.186716 0.445222 0.148689 0.147044 0.163881 0.357696 0.249228 0.120884 0.338756 0.23799 0.174653 0.250582 0.359913 0.191932 0.245742 0.165763 0.21248 0.261455 0.23596 0.279354 0.21639 0.253048 0.326013 0.274521 0.445528 0.214659 0.312161 0.156116 0.347797 0.169771 0.389098 0.294694 0.211097 0.264628 0.183466 0.236005 0.21618 0.347066 0.211115 0.320671 0.250266 0.159689 0.2296 0.207664 0.384994 0.242294 0.265087 0.151677 0.165181 0.309476 0.294808 0.205836 0.211114 0.158838 0.213561 0.259251 0.319854 0.186238 0.230296 0.389527 0.157746 0.135929 0.390339 0.190216 0.38471 0.270537 0.177239 0.233511 0.282228 0.313943 0.304133 0.387854 0.243608 0.322218 0.206424 0.203452 0.310101 0.262911 0.290483 0.191684 0.202841 0.20429 0.132531 0.192968 0.189274 0.198164 0.189636 0.292542 0.220207 0.224159 0.221258 0.249942 0.125013 0.172764 0.319641 0.190076 0.263651 0.310684 0.088757 0.240359 0.174358 0.198325 0.251556 0.280211 0.292003 0.18409 0.318648 0.137134 0.191295 0.16302 0.226551 0.24716 0.133549 0.258395 0.304624 0.264177 0.177385 0.274286 0.157228 0.299101 0.240278 0.194571 0.350429 0.212143 0.193681 0.266607 0.150228 0.204183 0.315672 0.193363 0.130019 0.202129 0.216326 0.213133 0.187089 0.190715 0.180934 0.183334 0.297964 0.277886 0.23907 0.207953 0.30527 0.176088 0.199663 0.289293 0.197831 0.258065 0.299011 0.179249 0.293239 0.204802 0.319032 0.211717 0.251688 0.242399 0.202651 0.198063 0.195439 0.343659 0.191462 0.302429 0.16561 0.271476 0.254929 0.185234 0.306123 0.325064 0.176234 0.367567 0.139686 0.147258 0.217225 0.15128 0.315976 0.265362 0.310024 0.292406 0.410225 0.286 0.274103 0.21506 0.265749 0.214389 0.142559 0.144568 0.117615 0.155886 0.255398 0.195037 0.192702 0.184238 0.224671 0.348663 0.150836 0.447537 0.163529 0.160546 0.177936 0.246547 0.231999 0.242593 0.254422 0.226703 0.271841 0.311305 0.173312 0.31315 0.2734 0.281242 0.247339 0.195463 0.248582 0.24879 0.198124 0.258069 0.246419 0.32682 0.234704 0.182361 0.242505 0.228687 0.18097 0.203595 0.186019 0.180037 0.387063 0.17064 0.300626 0.270758 0.233519 0.12908 0.279874 0.240328 0.201165 0.226879 0.17419 0.183829 0.183269 0.312356 0.235775 0.298285 0.175933 0.248909 0.317036 0.217625 0.287615 0.24738 0.281576 0.139421 0.402637 0.192828 0.240097 0.22426 0.300398 0.26835 0.203136 0.274225 0.281552 0.212764 0.223405 0.204213 0.225393 0.365084 0.284772 0.164882 0.167006 0.141252 0.334749 0.191537 0.347148 0.261467 0.278839 0.16289 0.156079 0.163327 0.453189 0.242531 0.211213 0.377327 0.129425 0.338062 0.178311 0.151441 0.177396 0.141849 0.260319 0.257188 0.253319 0.131483 0.17845 0.295193 0.291989 0.196744 0.245737 0.243616 0.238327 0.140141 0.299697 0.175678 0.222671 0.313991 0.229781 0.166001 0.267725 0.250766 0.178963 0.202857 0.265851 0.252263 0.186411 0.171526 0.203046 0.201325 0.219007 0.168171 0.150266 0.19056 0.225388 0.302085 0.253565 0.229089 0.169416 0.185946 0.237059 0.232741 0.207615 0.210586 0.271716 0.126831 0.143729 0.256381 0.153556 0.136276 0.346245 0.316635 0.197651 0.190807 0.180935 0.419781 0.29817 0.181433 0.383963 0.175703 0.290882 0.175076 0.267717 0.205832 0.21102 0.278328 0.253905 0.301712 0.207115 0.129098 0.230993 0.242539 0.420359 0.189212 0.230724 0.102015 0.185154 0.24309 0.200898 0.181461 0.36374 0.23083 0.242841 0.226207 0.13981 0.218217 0.323039 0.135569 0.356893 0.301319 0.204239 0.286343 0.214909 0.227409 0.291323 0.343346 0.194592 0.237034 0.260758 0.275149 0.188387 0.265764 0.289002 0.192824 0.264597 0.193837 0.107146 0.327921 0.239866 0.150128 0.162846 0.368589 0.224139 0.23766 0.148874 0.200754 0.37082 0.32886 0.23482 0.185807 0.263471 0.21654 0.294607 0.431925 0.141704 0.222322 0.254515 0.140725 0.243017 0.331751 0.225148 0.243863 0.277493 0.209899 0.14073 0.181801 0.158516 0.21906 0.143209 0.483481 0.184439 0.244044 0.138603 0.209734 0.183711 0.234662 0.206754 0.200536 0.37852 0.295702 0.153963 0.229749 0.37156 0.15438 0.318444 0.226383 0.174339 0.280018 0.298359 0.207417 0.191352 0.204486 0.25699 0.18336 0.203165 0.268912 0.124346 0.170587 0.208515 0.255785 0.189915 0.193742 0.265541 0.196593 0.44024 0.193236 0.210464 0.178401 0.267437 0.261669 0.151666 0.231747 0.429297 0.186123 0.247683 0.174011 0.20735 0.182367 0.206227 0.288376 0.215381 0.324549 0.293896 0.26096 0.168642 0.221237 0.138552 0.244411 0.242515 0.282381 0.200742 0.159237 0.173538 0.155255 0.286793 0.190673 0.212543 0.142557 0.163852 0.324514 0.226359 0.286034 0.265175 0.167512 0.1669 0.175965 0.21887 0.270963 0.213383 0.275588 0.344282 0.197601 0.23752 0.198386 0.22746 0.194644 0.270187 0.179387 0.291081 0.264052 0.313068 0.229684 0.220837 0.238326 0.235412 0.176756 0.174808 0.186154 0.236009 0.290265 0.223216 0.387859 0.131766 0.197939 0.319445 0.206208 0.392324 0.242403 0.118289 0.206046 0.187751 0.305221 0.250805 0.364331 0.267645 0.226781 0.159105 0.16499 0.26336 0.328654 0.346171 0.255116 0.160301 0.209395 0.281502 0.150833 0.336398 0.313434 0.329376 0.377569 0.20131 0.147744 0.262209 0.13104 0.223523 0.327143 0.20418 0.280188 0.291936 0.22639 0.212777 0.237577 0.240276 0.162601 0.329284 0.246391 0.257846 0.275798 0.16553 0.224596 0.214602 0.108482 0.122975 0.190303 0.1887 0.259566 0.217448 0.163838 0.402514 0.30647 0.231689 0.21686 0.368648 0.264375 0.150326 0.280215 0.33134 0.175181 0.204698 0.259924 0.222497 0.336779 0.184718 0.157597 0.210996 0.33826 0.285073 0.34046 0.293602 0.292515 0.27316 0.198009 0.236383 0.191815 0.289511 0.145669 0.149272 0.235195 0.274644 0.111119 0.377749 0.256948 0.319047 0.240941 0.127335 0.258226 0.343787 0.294881 0.196334 0.280911 0.190064 0.154285 0.249943 0.111489 0.181615 0.228158 0.211736 0.284702 0.197705 0.441521 0.195201 0.334151 0.268518 0.303251 0.306536 0.331138 0.184274 0.191848 0.372185 0.18715 0.178561 0.20491 0.283165 0.379047 0.169138 0.186049 0.191355 0.184844 0.170478 0.278305 0.198398 0.270559 0.147733 0.248628 0.22569 0.274181 0.268129 0.127316 0.232541 0.267888 0.190538 0.205132 0.196855 0.336186 0.26457 0.154453 0.314116 0.252485 0.200411 0.257389 0.181853 0.215055 0.372672 0.159009 0.460942 0.21763 0.3326 0.280303 0.248211 0.28108 0.260293 0.158741 0.20304 0.199559 0.264671 0.10865 0.324939 0.208012 0.348151 0.212992 0.2751 0.21011 0.177989 0.340371 0.226067 0.301772 0.267218 0.232069 0.269174 0.262403 0.351137 0.249238 0.184369 0.213975 0.278517 0.235683 0.159376 0.213882 0.31107 0.194552 0.415929 0.204497 0.22956 0.273111 0.129715 0.212743 0.345718 0.167485 0.441264 0.331441 0.243648 0.180703 0.321454 0.154117 0.24783 0.30002 0.266497 0.197979 0.16356 0.283438 0.205534 0.424048 0.245813 0.208954 0.275326 0.408611 0.249622 0.289372 0.231521 0.239632 0.201393 0.331421 0.169495 0.253773 0.126327 0.131185 0.157528 0.209044 0.135153 0.321253 0.22447 0.311219 0.440873 0.25541 0.297953 0.200286 0.233467 0.282013 0.151914 0.235911 0.24066 0.280032 0.194794 0.17807 0.28693 0.201501 0.258032 0.199711 0.238658 0.249272 0.113986 0.263756 0.160426 0.145255 0.226331 0.162874 0.151595 0.311545 0.22891 0.157676 0.372381 0.272108 0.228952 0.377745 0.371848 0.282495 0.366701 0.206004 0.167455 0.243538 0.127041 0.26622 0.194445 0.369589 0.127336 0.294265 0.23233 0.196285 0.178111 0.20804 0.207412 0.15939 0.164841 0.293223 0.239116 0.241927 0.182128 0.200257 0.392925 0.09585 0.277011 0.19302 0.168614 0.192435 0.267743 0.17074 0.168197 0.387367 0.15849 0.192768 0.29094 0.123732 0.253654 0.242846 0.200605 0.295268 0.238396 0.270247 0.277841 0.22903 0.286705 0.124793 0.166882 0.226336 0.320674 0.20695 0.268657 0.198469 0.205879 0.42506 0.157844 0.235845 0.258849 0.228159 0.237073 0.18473 0.168912 0.257687 0.22576 0.202789 0.296076 0.14461 0.210384 0.253176 0.213064 0.301611 0.136296 0.190611 0.129719 0.280651 0.175827 0.180506 0.163485 0.249197 0.298008 0.42439 0.203748 0.185139 0.194122 0.444812 0.217627 0.17685 0.152761 0.304344 0.273641 0.230422 0.199353 0.203255 0.252222 0.271944 0.25995 0.216406 0.212162 0.315117 0.206185 0.341405 0.265734 0.211672 0.281645 0.230907 0.181238 0.178727 0.353421 0.177786 0.293748 0.344134 0.267519 0.266972 0.356069 0.157132 0.249888 0.178505 0.241579 0.188815 0.240609 0.230228 0.271263 0.151587 0.172543 0.420536 0.219637 0.241604 0.139692 0.215672 0.31003 0.173712 0.405362 0.276985 0.24343 0.244759 0.220027 0.254031 0.227329 0.119022 0.18005 0.428474 0.209727 0.220906 0.201862 0.224455 0.176477 0.25072 0.195066 0.279537 0.191999 0.197355 0.187141 0.20208 0.151877 0.222726 0.231318 0.341245 0.218655 0.174157 0.20219 0.097281 0.21788 0.352508 0.090406 0.2596 0.147341 0.398229 0.167354 0.181078 0.162581 0.286929 0.266696 0.161263 0.165781 0.231949 0.285164 0.181632 0.331933 0.270385 0.426115 0.253774 0.196881 0.103544 0.210007 0.229234 0.200898 0.259641 0.172729 0.222599 0.337667 0.163576 0.289171 0.21922 0.232822 0.249138 0.161662 0.268214 0.166798 0.20127 0.179004 0.329356 0.183578 0.171552 0.247907 0.272317 0.29933 0.189544 0.244763 0.268099 0.31946 0.235941 0.177365 0.37809 0.18466 0.230496 0.33116 0.300799 0.206834 0.213464 0.223735 0.232158 0.228074 0.40477 0.211225 0.336018 0.203538 0.167296 0.142596 0.324765 0.200117 0.318541 0.293398 0.191155 0.181331 0.181682 0.196813 0.203476 0.320992 0.34201 0.27995 0.207427 0.181417 0.326642 0.258086 0.304588 0.325814 0.170579 0.271409 0.126369 0.253492 0.328468 0.19889 0.225717 0.09292 0.181577 0.193763 0.161723 0.175888 0.196816 0.322569 0.217792 0.22955 0.153886 0.273156 0.249384 0.312342 0.161403 0.23971 0.243869 0.224016 0.200884 0.267734 0.280596 0.166451 0.172362 0.237325 0.205192 0.200383 0.217297 0.186636 0.167022 0.19944 0.333785 0.172961 0.290578 0.153204 0.199777 0.172383 0.226133 0.349706 0.453844 0.301429 0.196569 0.288516 0.268282 0.25094 0.223975 0.192368 0.208871 0.147991 0.186904 0.378834 0.241599 0.242279 0.338148 0.202985 0.216276 0.206645 0.265583 0.261782 0.3235 0.222609 0.325371 0.276812 0.238677 0.193263 0.302357 0.384176 0.130408 0.267985 0.275394 0.239261 0.201292 0.192908 0.163195 0.268775 0.297527 0.386473 0.255482 0.143024 0.200622 0.20192 0.152857 0.334919 0.175805 0.191968 0.257786 0.242838 0.224994 0.168652 0.188705 0.372686 0.301237 0.214271 0.197083 0.273446 0.166654 0.312982 0.117791 0.160523 0.223314 0.344874 0.324368 0.260309 0.249657 0.255896 0.26608 0.365991 0.19582 0.2797 0.330746 0.164172 0.395057 0.123882 0.329342 0.191328 0.13972 0.128416 0.282948 0.318215 0.291334 0.33312 0.124429 0.179575 0.200245 0.188757 0.16641 0.322437 0.177082 0.243512 0.162092 0.276869 0.24671 0.240676 0.210606 0.346033 0.345078 0.340544 0.149935 0.222755 0.34059 0.235634 0.411319 0.272211 0.195498 0.159687 0.194439 0.328433 0.349274 0.435779 0.300383 0.159626 0.157511 0.40354 0.262588 0.161162 0.356404 0.258775 0.225231 0.134237 0.191709 0.21131 0.130036 0.150256 0.208149 0.302353 0.328351 0.407281 0.319536 0.217114 0.24602 0.159701 0.149808 0.222346 0.255441 0.213355 0.300626 0.225122 0.29984 0.251029 0.179887 0.227381 0.286893 0.25204 0.292665 0.232213 0.183487 0.254627 0.240664 0.273199 0.183619 0.195081 0.278324 0.214097 0.19816 0.155948 0.1499 0.215373 0.314024 0.180666 0.258613 0.318429 0.232011 0.31626 0.237138 0.221614 0.190927 0.301088 0.216936 0.099276 0.252303 0.140962 0.145021 0.238474 0.208927 0.268502 0.201336 0.309013 0.233675 0.257375 0.288891 0.207351 0.28951 0.194663 0.321175 0.250187 0.192218 0.259591 0.213861 0.131978 0.2554 0.218252 0.13025 0.342046 0.15638 0.353696 0.236078 0.157496 0.233963 0.221712 0.277247 0.307539 0.262752 0.254968 0.113966 0.360212 0.227047 0.205255 0.217922 0.253996 0.199294 0.244122 0.278394 0.273366 0.315034 0.261312 0.268346 0.173813 0.418683 0.156447 0.321702 0.216986 0.253102 0.243103 0.231738 0.226959 0.241027 0.139987 0.241648 0.200434 0.224784 0.246207 0.16624 0.161255 0.274403 0.260762 0.240069 0.381887 0.281215 0.299591 0.184737 0.200549 0.2276 0.245354 0.115874 0.2958 0.15587 0.252625 0.191411 0.240227 0.208048 0.220867 0.182651 0.314172 0.184599 0.199298 0.286831 0.193188 0.27514 0.223599 0.24109 0.224454 0.153253 0.354085 0.290956 0.166356 0.176684 0.152026 0.150564 0.146827 0.167858 0.152821 0.307615 0.230192 0.314625 0.177978 0.134134 0.342955 0.165484 0.187263 0.239363 0.208505 0.297543 0.291108 0.235031 0.166145 0.19086 0.291087 0.242302 0.268103 0.121746 0.221246 0.161116 0.226665 0.23058 0.272364 0.146539 0.292561 0.209006 0.299036 0.18716 0.178225 0.475067 0.192103 0.246041 0.222347 0.198257 0.277617 0.295317 0.330734 0.279825 0.155983 0.24814 0.347969 0.181631 0.191795 0.320126 0.277897 0.143625 0.229259 0.197924 0.216967 0.289048 0.145748 0.143255 0.184851 0.175606 0.323774 0.39006 0.120631 0.210932 0.291381 0.252915 0.217415 0.329352 0.313213 0.273236 0.129139 0.244101 0.199339 0.365694 0.187606 0.253345 0.251624 0.173361 0.268439 0.178831 0.26678 0.321886 0.24945 0.159504 0.27246 0.264804 0.195634 0.226736 0.195224 0.231929 0.170207 0.267698 0.217564 0.256474 0.256269 0.291163 0.29089 0.311739 0.250302 0.141006 0.339006 0.097808 0.29232 0.162476 0.210233 0.243343 0.228776 0.25933 0.164819 0.26286 0.280255 0.190462 0.273376 0.294625 0.214744 0.215956 0.229326 0.227174 0.246381 0.165428 0.206358 0.159576 0.18197 0.21729 0.261987 0.197658 0.187975 0.177736 0.152836 0.154232 0.241685 0.365471 0.161764 0.251372 0.270075 0.17172 0.316242 0.16298 0.17968 0.155713 0.303742 0.151052 0.241363 0.355575 0.24813 0.233237 0.277673 0.157859 0.206995 0.199851 0.362382 0.324692 0.280632 0.242855 0.294818 0.212852 0.137654 0.208245 0.349806 0.244284 0.298853 0.239382 0.338204 0.293436 0.231836 0.248504 0.248108 0.166593 0.270015 0.228858 0.382531 0.242716 0.328941 0.265578 0.209891 0.464243 0.279929 0.318466 0.318675 0.190722 0.138438 0.138866 0.290226 0.200209 0.192146 0.133888 0.268595 0.241417 0.181619 0.295298 0.147704 0.233148 0.381623 0.154685 0.200776 0.272329 0.185971 0.198999 0.223778 0.246835 0.314131 0.101513 0.232695 0.18026 0.218722 0.185553 0.251456 0.173059 0.16115 0.309133 0.205633 0.23382 0.192945 0.254208 0.215955 0.192462 0.273004 0.181087 0.19932 0.179532 0.23179 0.209626 0.174549 0.344108 0.188609 0.278161 0.229488 0.192631 0.176213 0.29581 0.388762 0.199572 0.141471 0.139969 0.326612 0.281038 0.169759 0.176582 0.204791 0.299035 0.268171 0.259304 0.210119 0.230358 0.259658 0.275492 0.298754 0.279697 0.200758 0.358125 0.241605 0.218593 0.267166 0.220701 0.238665 0.268 0.185778 0.244597 0.15268 0.333649 0.129635 0.195964 0.253282 0.241993 0.25664 0.34682 0.251184 0.396985 0.35925 0.20446 0.237751 0.2412 0.249095 0.277597 0.264785 0.317696 0.203243 0.146335 0.154609 0.201216 0.234114 0.186318 0.137144 0.315169 0.249675 0.257258 0.125148 0.29668 0.195354 0.229629 0.218084 0.359472 0.10946 0.207864 0.302714 0.153533 0.1861 0.262806 0.15891 0.296127 0.202217 0.181516 0.218736 0.247516 0.240319 0.207746 0.239314 0.328645 0.217816 0.200596 0.189782 0.19863 0.117772 0.352052 0.356886 0.268706 0.299535 0.301977 0.26863 0.186629 0.121302 0.326519 0.162098 0.29288 0.204406 0.119663 0.372314 0.27793 0.304912 0.146315 0.231425 0.23888 0.210119 0.344139 0.25159 0.195984 0.164867 0.234773 0.164782 0.161957 0.200635 0.311362 0.165114 0.138278 0.303915 0.16347 0.24749 0.200294 0.314936 0.269662 0.267422 0.232581 0.24667 0.251523 0.210364 0.13774 0.269239 0.352128 0.214665 0.22328 0.267379 0.212332 0.359702 0.211684 0.206213 0.324573 0.13891 0.150139 0.250431 0.190123 0.166131 0.20612 0.157311 0.301105 0.247052 0.217489 0.296792 0.176331 0.207701 0.151271 0.139976 0.170897 0.202592 0.158372 0.269387 0.313981 0.144082 0.393947 0.214453 0.239331 0.170201 0.213643 0.226999 0.29086 0.19445 0.164087 0.232344 0.234042 0.146208 0.372991 0.27716 0.184586 0.318902 0.336514 0.229271 0.12264 0.18242 0.364889 0.365021 0.245593 0.387806 0.281812 0.282463 0.275045 0.227098 0.379863 0.194685 0.273374 0.326785 0.273182 0.266637 0.271526 0.245007 0.302692 0.219787 0.242435 0.333099 0.262906 0.187084 0.177206 0.290464 0.214682 0.242584 0.233612 0.204756 0.288112 0.245238 0.248458 0.163262 0.248344 0.406899 0.17986 0.332851 0.230512 0.179981 0.202739 0.253908 0.163436 0.228924 0.204859 0.280152 0.313866 0.191614 0.459325 0.147027 0.19684 0.272991 0.23546 0.304755 0.196939 0.300616 0.16401 0.201323 0.222764 0.202436 0.335478 0.367052 0.222231 0.17619 0.216994 0.217893 0.243934 0.170129 0.241441 0.233152 0.238123 0.113562 0.32534 0.134679 0.221731 0.288803 0.211937 0.238672 0.279839 0.133376 0.22867 0.162877 0.188551 0.360954 0.252102 0.20929 0.178781 0.199974 0.205339 0.283435 0.203374 0.173516 0.366371 0.235728 0.261024 0.166097 0.245419 0.272608 0.279868 0.258661 0.230059 0.411808 0.165054 0.248477 0.304344 0.223246 0.197185 0.111705 0.202969 0.306996 0.250795 0.154076 0.309749 0.249814 0.196516 0.28005 0.280229 0.15957 0.193106 0.197056 0.189359 0.251477 0.132524 0.249214 0.238196 0.177097 0.308343 0.262402 0.169472 0.21217 0.28591 0.210981 0.228535 0.244216 0.182915 0.280609 0.120545 0.293697 0.162818 0.180662 0.196691 0.165971 0.248491 0.26977 0.298791 0.188067 0.257191 0.255899 0.311992 0.166251 0.242623 0.211095 0.220129 0.185352 0.364278 0.264943 0.28877 0.262754 0.301088 0.38903 0.198502 0.196052 0.289397 0.484756 0.222769 0.158654 0.19634 0.219344 0.351061 0.330042 0.231576 0.22116 0.208893 0.245017 0.148717 0.265133 0.22987 0.207085 0.180895 0.269457 0.141445 0.219303 0.125746 0.169098 0.248791 0.225552 0.278181 0.229112 0.2129 0.193465 0.197902 0.174548 0.336168 0.172827 0.283719 0.179986 0.441001 0.269635 0.179992 0.339356 0.204746 0.481891 0.433092 0.174169 0.206601 0.297956 0.244848 0.299297 0.171883 0.223646 0.204964 0.328324 0.150243 0.168011 0.190694 0.211146 0.327402 0.26133 0.198029 0.402721 0.32015 0.201762 0.282899 0.259025 0.181504 0.183495 0.330884 0.403594 0.304323 0.203356 0.242679 0.358643 0.365969 0.247855 0.278971 0.286218 0.140344 0.2326 0.177378 0.107505 0.22192 0.313195 0.249662 0.226558 0.329596 0.235021 0.261392 0.170945 0.221776 0.184411 0.297256 0.162801 0.179735 0.183539 0.262411 0.11764 0.226293 0.302204 0.206095 0.212405 0.272429 0.254948 0.250373 0.168959 0.287422 0.127831 0.241429 0.185522 0.170491 0.223114 0.315993 0.208152 0.135608 0.217047 0.210091 0.217443 0.267456 0.255399 0.133539 0.146609 0.148625 0.342143 0.153973 0.305663 0.248146 0.324247 0.480583 0.205371 0.352252 0.248132 0.271361 0.149666 0.166066 0.22776 0.157464 0.131343 0.200632 0.357463 0.23552 0.330112 0.143681 0.267284 0.142921 0.102561 0.202582 0.214182 0.260714 0.137265 0.21356 0.173791 0.276584 0.236861 0.214901 0.405456 0.28718 0.249494 0.215329 0.170559 0.184696 0.206938 0.202745 0.176986 0.280098 0.284193 0.256592 0.210717 0.176391 0.212355 0.179752 0.224631 0.20153 0.196059 0.259953 0.397948 0.166007 0.206848 0.29823 0.16137 0.37081 0.124626 0.140801 0.151351 0.457869 0.095296 0.407877 0.301253 0.344042 0.202801 0.224347 0.167033 0.171917 0.244425 0.368481 0.274259 0.16455 0.188338 0.166473 0.260002 0.227489 0.385067 0.238257 0.29848 0.160903 0.293283 0.223786 0.45449 0.312166 0.188499 0.213345 0.252709 0.241739 0.34157 0.251841 0.282935 0.247503 0.278192 0.228265 0.186759 0.166989 0.192469 0.215986 0.239813 0.177584 0.213252 0.15391 0.166439 0.274697 0.129117 0.170894 0.286307 0.287574 0.257685 0.144285 0.237463 0.186626 0.223071 0.46877 0.175756 0.186851 0.154076 0.235962 0.208293 0.138386 0.259382 0.150782 0.271781 0.10023 0.216003 0.293938 0.229279 0.39772 0.160784 0.185287 0.257577 0.273451 0.184477 0.195374 0.300646 0.224362 0.17679 0.288856 0.277477 0.277733 0.221874 0.288797 0.238764 0.229889 0.304939 0.318466 0.116485 0.242498 0.20455 0.329411 0.283387 0.213221 0.238907 0.233497 0.158246 0.337033 0.25968 0.306309 0.185304 0.134417 0.271497 0.229288 0.231784 0.289779 0.196335 0.139216 0.208856 0.3586 0.291266 0.309526 0.144637 0.217388 0.328197 0.267016 0.171029 0.247375 0.179604 0.205282 0.311256 0.18933 0.304491 0.146382 0.196348 0.269385 0.350421 0.236142 0.266633 0.300351 0.17756 0.32821 0.182441 0.257196 0.207876 0.285868 0.313672 0.265263 0.276883 0.149444 0.137043 0.153826 0.100785 0.182306 0.12298 0.225999 0.283299 0.229663 0.243546 0.180898 0.311683 0.243011 0.344524 0.239961 0.205737 0.253305 0.228106 0.257104 0.195947 0.201773 0.150244 0.197216 0.195369 0.238417 0.231611 0.25734 0.434599 0.216378 0.300062 0.17996 0.294788 0.130363 0.210653 0.191226 0.158813 0.139515 0.23573 0.293903 0.32386 0.265513 0.253324 0.189037 0.242871 0.291439 0.201673 0.320324 0.207218 0.097256 0.290256 0.178434 0.154896 0.203567 0.170963 0.209171 0.267402 0.334944 0.198926 0.197559 0.3539 0.174001 0.23701 0.136659 0.210017 0.227746 0.174649 0.327432 0.131812 0.173907 0.169121 0.245728 0.102962 0.244447 0.290412 0.261997 0.154073 0.262271 0.194868 0.277683 0.236513 0.290597 0.367656 0.184877 0.322065 0.294759 0.289353 0.142176 0.205517 0.188361 0.192143 0.243864 0.421146 0.203798 0.328863 0.297845 0.156264 0.296653 0.289526 0.150116 0.250196 0.227733 0.364857 0.181404 0.204845 0.258081 0.173435 0.217284 0.409493 0.19111 0.265598 0.201858 0.129704 0.298788 0.347585 0.286565 0.230514 0.170579 0.250203 0.329767 0.20729 0.198884 0.296673 0.146789 0.219464 0.177582 0.261113 0.284863 0.179575 0.218875 0.285113 0.267457 0.198401 0.300852 0.191097 0.248741 0.222807 0.282677 0.122128 0.16097 0.200511 0.310537 0.245171 0.261412 0.194629 0.368758 0.259514 0.377564 0.204435 0.232522 0.329717 0.259897 0.195991 0.17508 0.237948 0.194308 0.203112 0.151904 0.307271 0.14701 0.384102 0.20716 0.345286 0.347564 0.46987 0.210313 0.206393 0.304532 0.234605 0.249153 0.156441 0.262362 0.238647 0.175207 0.283825 0.222761 0.152844 0.315936 0.261983 0.133571 0.365974 0.2904 0.397867 0.327951 0.296269 0.244711 0.170972 0.176083 0.191568 0.282089 0.290983 0.249909 0.314096 0.212834 0.255752 0.22321 0.181205 0.171758 0.235227 0.319556 0.112468 0.206975 0.379276 0.136278 0.206278 0.242116 0.195468 0.206883 0.286594 0.15062 0.300638 0.145663 0.239141 0.192105 0.149972 0.310922 0.166146 0.208327 0.411285 0.205822 0.22639 0.301851 0.182321 0.3196 0.197652 0.20314 0.192185 0.146397 0.288301 0.182522 0.220736 0.279549 0.406222 0.255092 0.205004 0.23385 0.165358 0.210788 0.219974 0.211643 0.309108 0.189775 0.21495 0.224069 0.270028 0.177076 0.250817 0.285569 0.233097 0.138022 0.193281 0.372595 0.224983 0.325724 0.335975 0.196543 0.204522 0.169587 0.30201 0.212241 0.1982 0.283918 0.315446 0.172985 0.206613 0.368479 0.181768 0.245718 0.144167 0.183164 0.282982 0.246941 0.252668 0.181324 0.164732 0.320472 0.192042 0.17322 0.282579 0.387725 0.219637 0.160455 0.143459 0.272412 0.259103 0.331813 0.290442 0.332276 0.22538 0.216287 0.135249 0.275234 0.27031 0.20034 0.236472 0.183519 0.304182 0.158462 0.229578 0.296928 0.168836 0.296276 0.197477 0.192653 0.300017 0.27056 0.220644 0.240325 0.187911 0.193627 0.258758 0.201466 0.2867 0.250907 0.307535 0.39843 0.170255 0.268864 0.213671 0.19214 0.127572 0.360502 0.186808 0.197653 0.366121 0.269098 0.45957 0.274126 0.220922 0.256442 0.146625 0.133308 0.408723 0.190161 0.218703 0.251227 0.275377 0.286192 0.178095 0.203422 0.187316 0.233301 0.226181 0.270597 0.112449 0.155401 0.234215 0.137401 0.258329 0.354193 0.324927 0.204051 0.203498 0.197639 0.335312 0.276764 0.212208 0.271218 0.236945 0.34056 0.215678 0.318225 0.334039 0.202701 0.149282 0.240381 0.413342 0.259169 0.149769 0.2443 0.276154 0.138797 0.172381 0.169144 0.475419 0.260775 0.22135 0.256234 0.228378 0.447854 0.234824 0.335216 0.431328 0.337795 0.189116 0.143294 0.24748 0.181745 0.264245 0.253788 0.254618 0.262839 0.308601 0.245352 0.428514 0.162773 0.278665 0.215737 0.203829 0.154948 0.206548 0.170046 0.387916 0.186012 0.299902 0.199546 0.29728 0.185311 0.222407 0.220306 0.331661 0.21433 0.284835 0.256746 0.129631 0.244656 0.271317 0.146248 0.273873 0.212692 0.170569 0.198555 0.207712 0.233616 0.220087 0.162158 0.240926 0.235728 0.182928 0.183972 0.188596 0.120596 0.234323 0.224566 0.274815 0.425531 0.273419 0.174041 0.257608 0.218184 0.236128 0.284508 0.159289 0.259181 0.405347 0.204857 0.278933 0.185094 0.363881 0.39151 0.231009 0.150098 0.310971 0.140362 0.121999 0.268303 0.251871 0.195278 0.211197 0.266497 0.155323 0.163557 0.193417 0.280547 0.372581 0.17455 0.273472 0.299119 0.339743 0.152277 0.165014 0.318637 0.180154 0.258283 0.116509 0.237811 0.181517 0.29067 0.24614 0.204789 0.337648 0.208107 0.327982 0.32226 0.190913 0.243687 0.210922 0.285749 0.202149 0.328421 0.252231 0.327753 0.152756 0.214443 0.188685 0.204227 0.372078 0.205292 0.201795 0.345451 0.207361 0.157562 0.301224 0.263923 0.167015 0.174319 0.21736 0.232287 0.251791 0.200418 0.157811 0.254694 0.211928 0.303989 0.202148 0.099068 0.366078 0.176111 0.142178 0.328914 0.141849 0.250285 0.302979 0.217868 0.214635 0.213859 0.329322 0.159067 0.337456 0.273848 0.128203 0.270238 0.183189 0.253686 0.229772 0.177671 0.206346 0.184889 0.146278 0.287143 0.35102 0.289193 0.220052 0.168816 0.19598 0.373649 0.347907 0.283723 0.305959 0.261535 0.316815 0.210547 0.218157 0.443292 0.264067 0.254675 0.217622 0.238377 0.185526 0.239139 0.231843 0.208841 0.320618 0.199597 0.271804 0.265605 0.237988 0.194873 0.250607 0.274795 0.117624 0.216119 0.154311 0.237891 0.211423 0.254542 0.253161 0.267224 0.325878 0.249642 0.093623 0.183722 0.180116 0.078612 0.227954 0.177853 0.33283 0.232804 0.178788 0.253705 0.206639 0.290037 0.119572 0.197172 0.248893 0.416159 0.187377 0.230721 0.222301 0.138622 0.277851 0.310893 0.351784 0.292937 0.195367 0.298597 0.232037 0.262995 0.225355 0.23877 0.181258 0.252692 0.207008 0.169115 0.360351 0.225187 0.180811 0.352127 0.179914 0.175682 0.396229 0.320398 0.245831 0.161589 0.13709 0.232171 0.208331 0.193088 0.151956 0.172273 0.223349 0.16317 0.292758 0.257413 0.37407 0.208189 0.211409 0.247141 0.168323 0.236415 0.165761 0.164143 0.12812 0.186133 0.211862 0.213706 0.172492 0.243032 0.139258 0.162188 0.193454 0.193551 0.258954 0.211901 0.196544 0.153489 0.143035 0.201967 0.248471 0.202562 0.196304 0.288769 0.172753 0.218461 0.165723 0.176896 0.354043 0.281945 0.225216 0.12403 0.244431 0.240299 0.170681 0.355919 0.176348 0.178603 0.243298 0.303071 0.227474 0.20975 0.246631 0.239762 0.223106 0.356694 0.246979 0.258886 0.218313 0.345805 0.229701 0.318135 0.344295 0.187877 0.107119 0.232016 0.279558 0.218586 0.230938 0.194382 0.289172 0.23998 0.336662 0.145567 0.196316 0.176052 0.230562 0.212171 0.186896 0.179249 0.204573 0.223007 0.218733 0.1244 0.299261 0.231554 0.258841 0.22229 0.213725 0.186824 0.137976 0.150832 0.259466 0.314939 0.27137 0.179371 0.232718 0.140623 0.343293 0.167139 0.176719 0.248344 0.284055 0.242982 0.287049 0.310322 0.177478 0.242358 0.290533 0.121513 0.320112 0.264288 0.483725 0.234449 0.31395 0.174948 0.223445 0.241364 0.224523 0.230513 0.250933 0.151226 0.228623 0.226348 0.230199 0.165071 0.132516 0.320759 0.2365 0.168271 0.223998 0.194354 0.168902 0.327695 0.301529 0.29263 0.285953 0.288623 0.211794 0.257393 0.287105 0.274024 0.218661 0.145341 0.201168 0.216813 0.199562 0.292249 0.272213 0.200961 0.229191 0.317682 0.25617 0.246573 0.318992 0.274725 0.257768 0.214705 0.185266 0.426108 0.16897 0.332602 0.273543 0.145256 0.359014 0.189262 0.244355 0.299851 0.235575 0.229794 0.109743 0.189397 0.160897 0.207602 0.362125 0.185556 0.181118 0.278834 0.24467 0.351542 0.221864 0.140539 0.133917 0.382256 0.2202 0.316134 0.246721 0.183226 0.25874 0.358989 0.115081 0.276438 0.396209 0.25328 0.226343 0.247146 0.226953 0.209816 0.147201 0.233684 0.187428 0.231571 0.324575 0.303218 0.181092 0.203573 0.300293 0.116954 0.22447 0.308324 0.24925 0.196377 0.187606 0.314816 0.161609 0.198797 0.316276 0.243866 0.174466 0.276936 0.194829 0.317772 0.20299 0.143131 0.178104 0.141488 0.156807 0.263277 0.266711 0.213943 0.105625 0.267551 0.160165 0.151982 0.329209 0.217638 0.287799 0.292072 0.161295 0.288754 0.364118 0.164795 0.250314 0.183216 0.178878 0.162963 0.270569 0.192462 0.29505 0.168833 0.202001 0.321122 0.216598 0.220333 0.346158 0.315494 0.279422 0.162849 0.175327 0.149785 0.21359 0.232471 0.154803 0.257543 0.318759 0.161945 0.225865 0.242191 0.168412 0.209588 0.140042 0.209697 0.317882 0.261716 0.18008 0.215507 0.262738 0.218508 0.201799 0.190711 0.239999 0.219843 0.285453 0.144887 0.219508 0.307538 0.245397 0.1924 0.218096 0.214487 0.25883 0.129484 0.375384 0.262373 0.278978 0.206149 0.57343 0.211639 0.190419 0.245505 0.195818 0.262269 0.211356 0.277929 0.264895 0.187453 0.334669 0.28692 0.181052 0.248395 0.154359 0.196983 0.247032 0.158852 0.24399 0.164502 0.221208 0.175485 0.175297 0.243375 0.240346 0.163583 0.222927 0.188991 0.265056 0.221126 0.355919 0.298671 0.193412 0.335311 0.392459 0.211819 0.291024 0.409595 0.398341 0.238604 0.219714 0.241875 0.147113 0.280904 0.301698 0.362294 0.301867 0.132401 0.205738 0.20957 0.176448 0.273237 0.28294 0.168832 0.158425 0.262506 0.258073 0.336972 0.107351 0.174282 0.157798 0.206856 0.23436 0.232816 0.228282 0.222707 0.450354 0.155762 0.247049 0.110245 0.189385 0.147596 0.152817 0.126842 0.288697 0.216817 0.199641 0.192384 0.190316 0.2805 0.249709 0.148768 0.204618 0.106889 0.218001 0.318609 0.198682 0.312348 0.280001 0.411912 0.217235 0.219864 0.187718 0.311124 0.192153 0.230128 0.181255 0.308026 0.230804 0.407429 0.411445 0.267171 0.206949 0.132194 0.256136 0.278611 0.312986 0.228934 0.237059 0.192814 0.143646 0.215443 0.19303 0.142339 0.173268 0.169647 0.239064 0.264466 0.307347 0.157479 0.104411 0.210903 0.169387 0.24092 0.219217 0.155337 0.176125 0.180414 0.308152 0.226139 0.16228 0.208501 0.115518 0.318878 0.205041 0.298406 0.20187 0.259385 0.170035 0.253983 0.194686 0.224333 0.173548 0.160304 0.17543 0.206602 0.261433 0.377755 0.200489 0.205674 0.320363 0.252233 0.307394 0.278314 0.22312 0.265816 0.227126 0.176473 0.22854 0.281666 0.280052 0.209204 0.318277 0.305275 0.248618 0.31721 0.146481 0.223557 0.239997 0.243711 0.191419 0.204933 0.15187 0.13362 0.229373 0.193899 0.22785 0.245819 0.197019 0.227772 0.203059 0.180931 0.352873 0.19751 0.309126 0.211564 0.213527 0.157255 0.171918 0.318877 0.171493 0.207812 0.203397 0.25686 0.24004 0.177013 0.184365 0.152238 0.192503 0.186298 0.181875 0.153218 0.378931 0.180973 0.140513 0.162619 0.272607 0.267411 0.225671 0.332568 0.289326 0.164626 0.312384 0.212656 0.25544 0.179152 0.206202 0.235487 0.312598 0.113174 0.330416 0.255311 0.232443 0.179704 0.212185 0.139802 0.276429 0.283288 0.249785 0.369724 0.362965 0.260891 0.312728 0.353578 0.370444 0.222747 0.178686 0.315741 0.230377 0.474646 0.237244 0.267748 0.169699 0.135141 0.184835 0.201797 0.389685 0.354476 0.195883 0.230574 0.337432 0.224346 0.159599 0.26765 0.206691 0.149239 0.144302 0.224813 0.225883 0.14529 0.147883 0.198672 0.134104 0.253283 0.206287 0.09399 0.2083 0.20538 0.295263 0.262736 0.295073 0.162954 0.235518 0.233623 0.231481 0.245344 0.305358 0.223751 0.19038 0.248054 0.20119 0.255288 0.201072 0.188476 0.228838 0.214773 0.179184 0.377888 0.260224 0.269867 0.185782 0.293805 0.176608 0.200648 0.330029 0.129786 0.296138 0.149188 0.27812 0.275165 0.348778 0.203685 0.260736 0.285243 0.23684 0.298154 0.173044 0.211003 0.231158 0.195977 0.236936 0.182084 0.231485 0.323367 0.22388 0.241822 0.235921 0.113847 0.260671 0.290321 0.297017 0.230745 0.201571 0.30694 0.241266 0.215652 0.283819 0.143789 0.213754 0.389706 0.154725 0.23991 0.204904 0.273444 0.149949 0.135142 0.253675 0.225648 0.380227 0.192146 0.233102 0.184202 0.217825 0.236329 0.181627 0.454891 0.161927 0.1592 0.157666 0.347585 0.156051 0.209182 0.25128 0.166378 0.298468 0.194377 0.176571 0.211509 0.231705 0.264609 0.155302 0.172301 0.282449 0.176905 0.249533 0.201096 0.138529 0.288208 0.299349 0.210176 0.205119 0.432264 0.175083 0.190869 0.163151 0.281366 0.188402 0.171338 0.25819 0.191234 0.176454 0.260603 0.229701 0.356656 0.113973 0.201057 0.22803 0.18339 0.273918 0.227666 0.247315 0.314312 0.18407 0.200408 0.264336 0.276203 0.139534 0.294505 0.248409 0.350361 0.306904 0.212264 0.218187 0.23375 0.22048 0.260625 0.23084 0.199899 0.255271 0.192171 0.193543 0.178921 0.203456 0.283222 0.281878 0.201971 0.168182 0.275732 0.276384 0.226384 0.216236 0.208957 0.187665 0.304836 0.13869 0.336789 0.171226 0.249716 0.216671 0.253936 0.105585 0.180261 0.223877 0.216081 0.323202 0.262669 0.199581 0.167325 0.180204 0.373581 0.168741 0.281336 0.281879 0.232981 0.224643 0.249726 0.226327 0.31049 0.209165 0.283734 0.254965 0.209651 0.221879 0.282545 0.139746 0.340768 0.253561 0.294686 0.201648 0.285292 0.301138 0.147446 0.144627 0.286067 0.221545 0.215556 0.217961 0.105586 0.302931 0.203933 0.367366 0.179709 0.334305 0.313176 0.27887 0.263538 0.321511 0.145864 0.171689 0.165459 0.175551 0.227035 0.188832 0.232107 0.316467 0.233916 0.120698 0.363529 0.221922 0.221831 0.30065 0.209257 0.41412 0.317894 0.221456 0.273877 0.253773 0.165222 0.171276 0.218613 0.296617 0.260178 0.256571 0.187888 0.21551 0.286122 0.256674 0.216412 0.256223 0.245787 0.217718 0.14992 0.347175 0.179808 0.155555 0.213583 0.237474 0.19122 0.163968 0.139691 0.219942 0.23684 0.30821 0.244353 0.186579 0.219061 0.172175 0.261267 0.28928 0.187181 0.16936 0.189264 0.182858 0.274821 0.278254 0.206276 0.150325 0.246433 0.217603 0.369708 0.242243 0.149068 0.328911 0.338368 0.194754 0.252164 0.251873 0.328274 0.320543 0.177499 0.161674 0.228372 0.290234 0.205207 0.096317 0.096923 0.166774 0.198115 0.337292 0.411942 0.20185 0.415673 0.188339 0.253191 0.264021 0.179869 0.313744 0.203111 0.205856 0.187099 0.290727 0.186336 0.140581 0.161433 0.293508 0.239164 0.284227 0.19191 0.198167 0.346242 0.277354 0.317815 0.150125 0.150031 0.242235 0.2124 0.220849 0.22453 0.173648 0.185128 0.305056 0.342646 0.190039 0.20327 0.254343 0.176897 0.351528 0.187545 0.16939 0.290755 0.23238 0.386772 0.193188 0.243639 0.288521 0.165123 0.187589 0.268054 0.250645 0.226497 0.262107 0.234506 0.364979 0.240439 0.188383 0.240766 0.340306 0.122203 0.314658 0.211277 0.370842 0.240466 0.221485 0.286267 0.208991 0.293237 0.222548 0.16546 0.228715 0.263084 0.261818 0.259896 0.19362 0.218375 0.133858 0.148432 0.267206 0.221895 0.317284 0.195564 0.123355 0.283955 0.255138 0.110295 0.189698 0.332761 0.143846 0.287928 0.230891 0.287728 0.164279 0.198803 0.178988 0.14794 0.189008 0.211292 0.15966 0.188058 0.14718 0.267073 0.17691 0.325716 0.146831 0.24635 0.187832 0.25827 0.240705 0.210751 0.262487 0.225416 0.186345 0.252703 0.134982 0.300825 0.276575 0.166599 0.198953 0.25383 0.144743 0.168534 0.220374 0.337937 0.25602 0.179536 0.262971 0.194051 0.286624 0.249658 0.392265 0.161107 0.156855 0.160323 0.265675 0.177595 0.147648 0.289967 0.218854 0.223721 0.294621 0.271111 0.235376 0.260131 0.124378 0.184298 0.160937 0.217599 0.233755 0.222669 0.164931 0.281534 0.194123 0.235025 0.236608 0.22311 0.164263 0.363247 0.197518 0.163303 0.281008 0.135076 0.14275 0.222687 0.196614 0.267278 0.234485 0.242206 0.252618 0.274046 0.218027 0.206107 0.312394 0.277794 0.22531 0.163027 0.331618 0.293951 0.184228 0.128735 0.348641 0.216378 0.1898 0.20313 0.180804 0.158894 0.320081 0.364895 0.3629 0.2487 0.279042 0.260491 0.162024 0.237883 0.196451 0.188068 0.179863 0.237542 0.360855 0.238377 0.19153 0.122208 0.137615 0.204826 0.275072 0.225806 0.17256 0.173007 0.256453 0.151667 0.324903 0.31664 0.128612 0.196453 0.234443 0.199965 0.394547 0.212845 0.265476 0.233559 0.211794 0.29441 0.230533 0.184569 0.265988 0.241097 0.284992 0.237238 0.262641 0.229829 0.291843 0.157822 0.350544 0.218213 0.210093 0.176871 0.258043 0.247164 0.247373 0.148718 0.248885 0.169126 0.369148 0.265683 0.280829 0.18439 0.221549 0.263197 0.168277 0.2243 0.204985 0.348101 0.197224 0.286062 0.18196 0.266567 0.275779 0.142841 0.191102 0.14934 0.259451 0.271102 0.117269 0.27038 0.29433 0.164367 0.285016 0.273348 0.188844 0.215176 0.25554 0.229048 0.280508 0.191844 0.436136 0.230751 0.20519 0.233097 0.197873 0.260932 0.168946 0.160776 0.242441 0.200455 0.152386 0.307787 0.171557 0.267561 0.29757 0.305937 0.289017 0.291138 0.198067 0.322871 0.263379 0.309677 0.330828 0.245583 0.283254 0.207052 0.365913 0.328877 0.130365 0.18811 0.174673 0.165524 0.338779 0.303503 0.317472 0.280905 0.168499 0.252201 0.285991 0.25531 0.216517 0.233313 0.309434 0.233226 0.195712 0.210084 0.325172 0.237565 0.165392 0.305682 0.202917 0.208235 0.24329 0.322361 0.135921 0.181967 0.230965 0.360481 0.293698 0.24178 0.271484 0.320058 0.34749 0.23694 0.385248 0.233564 0.297573 0.317518 0.226677 0.27316 0.252444 0.339385 0.184736 0.195619 0.251585 0.236906 0.233884 0.267589 0.20817 0.170079 0.281496 0.120064 0.329812 0.206924 0.146162 0.205298 0.217038 0.234019 0.254633 0.227062 0.245388 0.328424 0.311479 0.142176 0.181628 0.138289 0.158991 0.262102 0.12328 0.189979 0.22222 0.261191 0.166804 0.253332 0.170004 0.248274 0.393029 0.181013 0.216186 0.295772 0.253775 0.203086 0.265295 0.384522 0.184771 0.193095 0.32504 0.276179 0.338941 0.234012 0.183634 0.175166 0.260932 0.399255 0.256443 0.311584 0.190017 0.155244 0.320888 0.323514 0.20086 0.186042 0.339987 0.192095 0.258159 0.213614 0.202917 0.310802 0.203849 0.197933 0.183769 0.233415 0.2151 0.252444 0.392226 0.205146 0.285835 0.197835 0.229252 0.265241 0.169492 0.286018 0.238581 0.175916 0.227624 0.226692 0.227842 0.208836 0.17024 0.209223 0.22123 0.175364 0.185952 0.25544 0.275543 0.307232 0.243666 0.17951 0.178917 0.183105 0.198419 0.233271 0.29102 0.178057 0.222383 0.208972 0.139766 0.326552 0.24557 0.194348 0.297482 0.172752 0.295478 0.219948 0.242703 0.213697 0.252859 0.271134 0.254819 0.208141 0.362521 0.137284 0.243985 0.160957 0.222265 0.218272 0.211254 0.245384 0.209079 0.327287 0.302594 0.258271 0.383944 0.139286 0.140028 0.277401 0.211354 0.249555 0.365651 0.337383 0.180009 0.290767 0.160593 0.224411 0.193801 0.324467 0.199166 0.402959 0.162498 0.175395 0.181263 0.395091 0.220226 0.169314 0.175683 0.177639 0.211039 0.295151 0.35944 0.21548 0.243372 0.361026 0.335721 0.355029 0.159981 0.214816 0.180336 0.246983 0.315208 0.185906 0.150042 0.451832 0.221587 0.107574 0.305065 0.293773 0.202281 0.194372 0.131766 0.22024 0.388331 0.226732 0.257727 0.165111 0.222774 0.235728 0.21615 0.327074 0.191441 0.313468 0.204664 0.136088 0.225606 0.221505 0.278928 0.20319 0.231705 0.180929 0.326466 0.40486 0.144902 0.189965 0.22916 0.263785 0.265095 0.22828 0.157363 0.223822 0.214003 0.203318 0.237204 0.155265 0.30278 0.1188 0.232516 0.294035 0.228759 0.229183 0.157195 0.197249 0.231037 0.295182 0.272821 0.325958 0.276465 0.188995 0.227841 0.23922 0.152097 0.241856 0.298951 0.2934 0.246475 0.269783 0.261763 0.257785 0.268881 0.107015 0.288193 0.178483 0.175907 0.242117 0.333046 0.380888 0.152406 0.173024 0.26808 0.271055 0.191567 0.213477 0.237378 0.230065 0.321374 0.323679 0.213069 0.230949 0.236823 0.323459 0.244842 0.353522 0.21205 0.192949 0.228411 0.227057 0.25422 0.354528 0.201993 0.299728 0.155009 0.312412 0.292328 0.275059 0.160173 0.188964 0.202159 0.256585 0.3633 0.344512 0.262084 0.236219 0.237747 0.306274 0.212228 0.438358 0.112795 0.247097 0.151439 0.21509 0.198316 0.296031 0.277176 0.227099 0.145808 0.201568 0.190658 0.219097 0.324417 0.30859 0.318367 0.219486 0.167673 0.291808 0.215118 0.158739 0.209449 0.157263 0.181464 0.312884 0.119797 0.213638 0.312258 0.252694 0.265557 0.217709 0.18953 0.248465 0.338782 0.261488 0.319132 0.227642 0.312632 0.241967 0.254745 0.284749 0.267524 0.237731 0.209747 0.229352 0.224634 0.277762 0.228036 0.258337 0.141365 0.217024 0.206375 0.23793 0.158722 0.172782 0.286496 0.196462 0.168714 0.169693 0.152005 0.158303 0.193959 0.173069 0.249664 0.233551 0.22768 0.2536 0.211877 0.253035 0.339658 0.193909 0.131782 0.258821 0.125006 0.286342 0.173783 0.237174 0.176616 0.268311 0.370543 0.26529 0.337803 0.138543 0.226629 0.221356 0.233843 0.251216 0.214658 0.227195 0.273457 0.145681 0.240214 0.128754 0.197965 0.234361 0.219166 0.223671 0.323232 0.287072 0.293149 0.244175 0.240139 0.375695 0.281618 0.301882 0.229626 0.240842 0.179501 0.159354 0.213158 0.160264 0.166697 0.153247 0.23465 0.222801 0.361684 0.202523 0.196708 0.217693 0.132122 0.190118 0.308054 0.342828 0.166475 0.344271 0.190081 0.258594 0.283716 0.202395 0.181088 0.182929 0.099272 0.088132 0.23066 0.258531 0.230062 0.174034 0.354983 0.158098 0.191794 0.288545 0.153644 0.232533 0.266652 0.257288 0.245977 0.202735 0.375901 0.153577 0.20906 0.284284 0.339942 0.207434 0.299789 0.154792 0.21453 0.297585 0.172481 0.305305 0.382738 0.224605 0.301592 0.399922 0.286546 0.247501 0.274322 0.242236 0.219618 0.360207 0.318057 0.152466 0.126845 0.177822 0.206594 0.202086 0.331322 0.240933 0.546108 0.336409 0.178582 0.231941 0.235575 0.215407 0.300974 0.187724 0.272827 0.228995 0.239921 0.182915 0.224983 0.256291 0.341513 0.168362 0.333244 0.181707 0.204166 0.175424 0.33274 0.388901 0.170898 0.234153 0.349117 0.14585 0.323969 0.316939 0.363912 0.231359 0.223965 0.1571 0.328777 0.188446 0.42519 0.200688 0.187142 0.205781 0.463341 0.327279 0.239034 0.195131 0.309144 0.167438 0.244939 0.176582 0.252849 0.151425 0.189843 0.175391 0.1545 0.293925 0.202239 0.398887 0.198551 0.163039 0.276384 0.235838 0.153076 0.155046 0.352329 0.223524 0.20443 0.213073 0.207427 0.2533 0.142344 0.27262 0.164158 0.232146 0.26401 0.253533 0.182672 0.16812 0.158855 0.339989 0.314212 0.250583 0.236343 0.247241 0.263504 0.310648 0.345063 0.189617 0.207623 0.175922 0.141294 0.182237 0.224623 0.185913 0.137854 0.115962 0.281536 0.180581 0.260174 0.196171 0.155052 0.182704 0.127095 0.418121 0.319038 0.224309 0.289752 0.215073 0.221506 0.232391 0.263044 0.165416 0.302648 0.247691 0.200958 0.265989 0.136426 0.340036 0.158261 0.337673 0.227929 0.227461 0.184889 0.351761 0.157257 0.242077 0.196385 0.223321 0.280406 0.214515 0.116964 0.232726 0.31043 0.191885 0.195933 0.253454 0.274834 0.340076 0.318512 0.213972 0.245841 0.310296 0.293394 0.231538 0.252704 0.292625 0.181791 0.298268 0.263822 0.232157 0.199623 0.364757 0.159217 0.225799 0.184406 0.217883 0.274138 0.177304 0.238262 0.283521 0.308749 0.346083 0.364581 0.208579 0.286943 0.173928 0.280191 0.141788 0.186109 0.171277 0.260378 0.205556 0.255064 0.3067 0.274089 0.172942 0.19878 0.131784 0.280486 0.147263 0.18707 0.183521 0.377523 0.19217 0.407848 0.214253 0.14025 0.323552 0.297242 0.216875 0.17997 0.243009 0.259801 0.309077 0.202857 0.145277 0.234963 0.24687 0.159745 0.208554 0.172612 0.215327 0.182444 0.456791 0.222702 0.235358 0.319443 0.356413 0.178515 0.233729 0.205545 0.342652 0.272036 0.23499 0.233874 0.186421 0.21981 0.241659 0.108387 0.174801 0.152115 0.188412 0.200434 0.174443 0.221383 0.275298 0.210181 0.152903 0.193195 0.413373 0.30148 0.229603 0.324866 0.210424 0.286354 0.137611 0.376203 0.252775 0.280486 0.304114 0.223486 0.274048 0.292001 0.245615 0.175476 0.431564 0.269682 0.202141 0.126784 0.274482 0.32128 0.156445 0.204222 0.383347 0.279196 0.356414 0.236676 0.302855 0.214162 0.326177 0.228256 0.208359 0.17542 0.152304 0.239845 0.266564 0.216028 0.197909 0.213693 0.287269 0.492666 0.307187 0.20614 0.216392 0.177683 0.253318 0.128926 0.220353 0.154543 0.165102 0.14449 0.171558 0.274295 0.156355 0.171097 0.276928 0.25023 0.168986 0.341203 0.287412 0.145344 0.248043 0.262807 0.225642 0.28677 0.232663 0.213311 0.221527 0.255948 0.207985 0.15721 0.249466 0.283267 0.296046 0.183373 0.239341 0.354833 0.321075 0.206462 0.15923 0.265128 0.336758 0.159419 0.181292 0.243611 0.25523 0.232039 0.208077 0.299989 0.278017 0.26253 0.169694 0.323809 0.320306 0.20805 0.1897 0.209147 0.263607 0.198953 0.30109 0.200048 0.282966 0.124897 0.207129 0.170203 0.243608 0.159068 0.134055 0.172694 0.284631 0.229964 0.201212 0.157472 0.143996 0.278326 0.241097 0.178358 0.248715 0.146239 0.15566 0.215223 0.164016 0.282377 0.32342 0.260367 0.161108 0.318592 0.277967 0.359882 0.120683 0.365965 0.260171 0.16833 0.399633 0.159025 0.372548 0.291663 0.423278 0.24818 0.189827 0.339972 0.136522 0.204198 0.27821 0.179318 0.181779 0.222238 0.313069 0.164834 0.188028 0.170161 0.239004 0.291541 0.312176 0.234958 0.291724 0.147143 0.221671 0.15793 0.139342 0.173013 0.161603 0.233039 0.259129 0.24847 0.216928 0.217868 0.165336 0.329212 0.363838 0.208998 0.258603 0.222602 0.144238 0.191664 0.22981 0.220616 0.25741 0.283971 0.186575 0.244489 0.136003 0.159323 0.187311 0.338528 0.264698 0.28664 0.224963 0.251627 0.242768 0.251366 0.16425 0.2579 0.201571 0.206997 0.25901 0.202013 0.179658 0.320974 0.366496 0.224227 0.17335 0.29656 0.236587 0.15534 0.205755 0.249384 0.345217 0.173942 0.205003 0.256466 0.233619 0.217485 0.212834 0.199363 0.208862 0.192502 0.171178 0.339945 0.29044 0.384432 0.195247 0.251139 0.22537 0.116694 0.146042 0.177142 0.175738 0.148918 0.237333 0.315024 0.163719 0.248844 0.188595 0.105708 0.278456 0.15784 0.298575 0.335053 0.179228 0.213692 0.231035 0.204862 0.144225 0.254947 0.222734 0.276816 0.199339 0.365204 0.177758 0.164681 0.158208 0.353513 0.320917 0.242317 0.246571 0.291497 0.192515 0.233395 0.219679 0.272306 0.192021 0.208864 0.134487 0.147858 0.136303 0.164506 0.246334 0.282239 0.179153 0.202847 0.175621 0.277718 0.219249 0.266052 0.197413 0.378211 0.176256 0.220042 0.172841 0.133124 0.270096 0.260561 0.145778 0.165559 0.292049 0.223743 0.220922 0.225789 0.173209 0.250425 0.19545 0.293944 0.278722 0.419454 0.174254 0.260942 0.100311 0.252104 0.265101 0.267808 0.19086 0.164659 0.319995 0.170478 0.160015 0.225594 0.240894 0.275018 0.23463 0.228662 0.164835 0.242045 0.188885 0.425381 0.240803 0.196594 0.368781 0.26071 0.282399 0.204763 0.253864 0.345782 0.175379 0.286707 0.14977 0.210763 0.254643 0.357566 0.19396 0.42428 0.293878 0.268496 0.221308 0.301642 0.234572 0.146122 0.227636 0.181061 0.245843 0.244883 0.16653 0.333005 0.147438 0.186114 0.314229 0.212377 0.240892 0.195982 0.28045 0.182621 0.167007 0.332829 0.279959 0.250907 0.237833 0.178629 0.197342 0.140274 0.117349 0.221629 0.190439 0.254998 0.213631 0.242695 0.270037 0.326446 0.141915 0.161909 0.216161 0.267268 0.22752 0.198617 0.190789 0.384093 0.332599 0.313002 0.187186 0.275668 0.188955 0.223866 0.293702 0.174935 0.181954 0.25359 0.244571 0.183375 0.357141 0.242502 0.253475 0.195558 0.274425 0.149907 0.165809 0.120554 0.166054 0.258781 0.151004 0.182312 0.253136 0.281942 0.324916 0.212256 0.161177 0.139227 0.182349 0.221994 0.214041 0.151497 0.137191 0.373582 0.207948 0.320691 0.226463 0.306149 0.385824 0.246024 0.243635 0.249359 0.150324 0.280549 0.168614 0.166259 0.267516 0.136624 0.146012 0.373485 0.33104 0.265974 0.254496 0.15526 0.173095 0.166597 0.180675 0.322075 0.383321 0.287371 0.256354 0.158914 0.183533 0.228564 0.226043 0.206837 0.122916 0.269877 0.193065 0.255538 0.269283 0.237061 0.320712 0.276288 0.457532 0.344373 0.23707 0.234642 0.334729 0.261416 0.280006 0.362218 0.314583 0.260686 0.23115 0.216959 0.229362 0.218675 0.22236 0.289637 0.215584 0.136622 0.164708 0.278397 0.139068 0.158521 0.376202 0.32436 0.255846 0.245119 0.200906 0.240273 0.250452 0.310471 0.168364 0.125644 0.197924 0.26715 0.188961 0.208738 0.301435 0.361899 0.23125 0.203284 0.241469 0.186017 0.310902 0.195992 0.15028 0.206231 0.301909 0.181142 0.16756 0.188861 0.178113 0.185195 0.19705 0.197056 0.191723 0.300292 0.263111 0.337332 0.213625 0.191513 0.161575 0.154749 0.314291 0.278393 0.235977 0.337481 0.178186 0.329933 0.261725 0.273471 0.182185 0.3576 0.333137 0.140827 0.246463 0.165003 0.238944 0.215906 0.196973 0.442972 0.221436 0.209307 0.225779 0.209395 0.241293 0.212727 0.330728 0.266532 0.264173 0.204205 0.182916 0.334667 0.325812 0.223101 0.21387 0.183812 0.239158 0.268947 0.216178 0.270881 0.346016 0.195879 0.240575 0.279173 0.307897 0.1867 0.314513 0.275511 0.215392 0.167248 0.229183 0.150276 0.196582 0.19698 0.309329 0.319065 0.17914 0.231746 0.164167 0.230287 0.217095 0.201973 0.193539 0.173235 0.266856 0.206535 0.260168 0.22954 0.171726 0.169431 0.166559 0.201029 0.249604 0.242089 0.405912 0.217494 0.180087 0.133497 0.227408 0.167909 0.240302 0.283219 0.258516 0.314547 0.192814 0.24813 0.115779 0.242864 0.224803 0.172579 0.212556 0.264853 0.313622 0.218863 0.223986 0.217268 0.15703 0.229944 0.290575 0.183852 0.261289 0.22344 0.134675 0.140652 0.22875 0.227235 0.168276 0.145793 0.262917 0.242001 0.37979 0.120969 0.22584 0.16868 0.146451 0.319126 0.215556 0.218881 0.239521 0.182446 0.259127 0.257163 0.186542 0.187199 0.193249 0.272311 0.202044 0.222369 0.219321 0.15177 0.28202 0.224091 0.291998 0.128795 0.206568 0.300192 0.242974 0.129432 0.168398 0.200226 0.226934 0.117876 0.171123 0.232549 0.351193 0.212998 0.310581 0.272461 0.262568 0.259847 0.333643 0.34778 0.266221 0.283277 0.284243 0.214946 0.528641 0.405368 0.239281 0.260244 0.23587 0.290594 0.197858 0.189754 0.352848 0.22951 0.25877 0.235025 0.246194 0.113373 0.177091 0.327042 0.143779 0.332569 0.253613 0.225049 0.245209 0.24339 0.275621 0.1924 0.190787 0.31351 0.18042 0.220124 0.322253 0.421941 0.169662 0.200259 0.13408 0.166461 0.150407 0.157273 0.307495 0.25796 0.31536 0.237969 0.301578 0.260971 0.235826 0.251548 0.136286 0.188777 0.209495 0.383601 0.178619 0.145433 0.171506 0.263104 0.259449 0.330113 0.221979 0.25345 0.347768 0.250107 0.139012 0.197373 0.163009 0.157959 0.447449 0.297804 0.110098 0.236904 0.237724 0.274078 0.217723 0.24598 0.319714 0.407719 0.251613 0.311596 0.366897 0.267304 0.236852 0.208653 0.211364 0.287659 0.263964 0.215193 0.180458 0.357823 0.296985 0.219213 0.207353 0.21737 0.183422 0.372372 0.190133 0.284544 0.182738 0.279242 0.264144 0.217859 0.308196 0.16855 0.211915 0.235962 0.360914 0.138721 0.225492 0.258684 0.302873 0.368894 0.186775 0.129187 0.260087 0.306775 0.25516 0.24179 0.216952 0.339488 0.274943 0.206807 0.422173 0.35049 0.192518 0.302805 0.159652 0.157953 0.345327 0.241538 0.4026 0.260296 0.190616 0.175279 0.247912 0.211426 0.179567 0.185163 0.23574 0.30108 0.230948 0.213497 0.291655 0.125145 0.212003 0.24689 0.15884 0.185338 0.270422 0.181587 0.286041 0.298866 0.184164 0.28583 0.318438 0.203712 0.256871 0.282089 0.280796 0.265989 0.18014 0.283098 0.262281 0.168354 0.286258 0.172304 0.393543 0.287219 0.213204 0.157688 0.254762 0.407985 0.254327 0.196957 0.188152 0.459075 0.158482 0.220608 0.342423 0.252881 0.309563 0.334738 0.240666 0.290904 0.244248 0.102516 0.192854 0.193303 0.18398 0.231311 0.184885 0.319139 0.188903 0.210831 0.145037 0.251494 0.187725 0.289866 0.294664 0.203354 0.365267 0.182759 0.296762 0.206088 0.332457 0.157248 0.233846 0.393769 0.196761 0.16087 0.271166 0.342262 0.195467 0.360683 0.178737 0.273618 0.125042 0.120738 0.273606 0.291985 0.183408 0.143146 0.235931 0.340994 0.201391 0.299101 0.11615 0.241726 0.22364 0.292379 0.232726 0.180041 0.26355 0.189008 0.211937 0.156671 0.33745 0.219405 0.189689 0.188423 0.232636 0.256912 0.248689 0.148485 0.234558 0.261409 0.294125 0.165677 0.223883 0.21845 0.286242 0.286686 0.261119 0.359385 0.123164 0.225418 0.16261 0.266507 0.282724 0.262261 0.11718 0.272593 0.169725 0.151553 0.218099 0.285996 0.171492 0.206453 0.153569 0.215611 0.28474 0.337075 0.228047 0.279195 0.233895 0.202152 0.252839 0.179486 0.328509 0.200028 0.184494 0.255197 0.353579 0.226864 0.244514 0.145101 0.286972 0.22275 0.282005 0.22316 0.253923 0.252804 0.314307 0.246982 0.221924 0.130803 0.25218 0.218189 0.242011 0.208654 0.248802 0.326837 0.280522 0.295164 0.133482 0.290623 0.192041 0.190284 0.491872 0.32763 0.211625 0.324548 0.286352 0.315107 0.279416 0.176044 0.198929 0.269939 0.196606 0.201491 0.251681 0.219921 0.336223 0.38202 0.227915 0.180232 0.169917 0.099272 0.413877 0.200852 0.18159 0.230233 0.351524 0.195126 0.172104 0.268887 0.129778 0.395065 0.167352 0.208375 0.194191 0.256811 0.334807 0.226568 0.343276 0.17038 0.253858 0.250749 0.234617 0.147359 0.244141 0.130968 0.182012 0.258505 0.198696 0.283615 0.200021 0.222001 0.280275 0.18579 0.207292 0.172255 0.158696 0.20947 0.316418 0.245602 0.19983 0.227355 0.263367 0.29919 0.286241 0.256039 0.256669 0.142914 0.203257 0.251206 0.237824 0.319048 0.161851 0.206522 0.118867 0.365379 0.280716 0.271485 0.303408 0.427703 0.368188 0.221742 0.192097 0.194193 0.138422 0.154473 0.2097 0.248737 0.297382 0.331101 0.189792 0.25822 0.212174 0.170473 0.217546 0.396493 0.382766 0.247712 0.221878 0.384057 0.31987 0.265083 0.255956 0.184275 0.311883 0.3143 0.173233 0.159104 0.207954 0.195965 0.267353 0.241271 0.161676 0.324755 0.173363 0.216503 0.229934 0.244216 0.183025 0.169782 0.186986 0.327616 0.19865 0.339035 0.362127 0.209655 0.270092 0.225771 0.304957 0.232159 0.322492 0.160367 0.21129 0.228736 0.301187 0.234348 0.146062 0.234711 0.253244 0.218837 0.187932 0.236982 0.199852 0.319351 0.227549 0.167438 0.36206 0.225741 0.154082 0.228015 0.148079 0.149186 0.263675 0.257662 0.210724 0.26023 0.189848 0.248419 0.290834 0.306647 0.203017 0.140074 0.214009 0.31097 0.322496 0.273376 0.167071 0.300443 0.20277 0.236063 0.268678 0.203071 0.207626 0.183191 0.380737 0.190255 0.187034 0.308961 0.192761 0.264321 0.315835 0.233369 0.191222 0.240981 0.151459 0.356858 0.25246 0.27716 0.224752 0.234594 0.252946 0.268078 0.269419 0.262404 0.28639 0.274473 0.246346 0.181287 0.261896 0.163237 0.178883 0.186498 0.273074 0.18267 0.250661 0.250203 0.21014 0.161608 0.249704 0.22072 0.146656 0.369776 0.271477 0.307741 0.15551 0.274653 0.313466 0.285069 0.214909 0.184095 0.27444 0.154252 0.267729 0.240622 0.27394 0.206135 0.156652 0.251868 0.350212 0.185711 0.209923 0.219993 0.185566 0.177247 0.243399 0.209985 0.220443 0.351177 0.140513 0.188501 0.202287 0.192101 0.263109 0.297 0.256032 0.265544 0.1216 0.204077 0.220117 0.333757 0.275762 0.124531 0.253478 0.313128 0.251992 0.234328 0.176548 0.277339 0.227404 0.306006 0.207021 0.168155 0.204717 0.259502 0.232773 0.29536 0.261761 0.391866 0.206795 0.201518 0.295343 0.291341 0.163518 0.132881 0.235925 0.322499 0.292043 0.215323 0.152311 0.251626 0.243687 0.149 0.217688 0.340785 0.120716 0.390645 0.263861 0.236478 0.33078 0.204547 0.25873 0.361964 0.299804 0.223005 0.162038 0.271333 0.148422
Exponential_Weibull 0.230347 0.308228 0.204545 0.331749 0.228132 0.220478 0.273842 0.448235 0.172393 0.161881 0.376545 0.188924 0.289196 0.301164 0.334944 0.206075 0.368479 0.358583 0.30587 0.213286 0.396336 0.264508 0.2323 0.234527 0.269605 0.211686 0.222507 0.235289 0.345617 0.298119 0.336142 0.248639 0.256737 0.252308 0.279153 0.185191 0.353858 0.168317 0.2273 0.065823 0.311691 0.24726 0.345622 0.211932 0.36428 0.290131 0.215653 0.345195 0.244121 0.319582 0.267831 0.334238 0.31698 0.199182 0.319463 0.343102 0.125034 0.235652 0.275513 0.265254 0.213532 0.12318 0.242485 0.245237 0.292937 0.204117 0.258896 0.305893 0.29234 0.31882 0.183769 0.243189 0.217731 0.32018 0.286999 0.333408 0.340482 0.264785 0.198065 0.283119 0.223444 0.136172 0.195501 0.201825 0.210003 0.294132 0.219265 0.22541 0.213938 0.151793 0.305107 0.273592 0.338859 0.262155 0.293538 0.130221 0.116886 0.231489 0.291424 0.203177 0.20762 0.295461 0.192489 0.303341 0.273077 0.298298 0.32993 0.280698 0.306284 0.288206 0.328454 0.256833 0.316285 0.335849 0.331204 0.295812 0.269866 0.146342 0.189217 0.261828 0.283773 0.366307 0.268144 0.305618 0.283073 0.176174 0.100005 0.346218 0.229932 0.168976 0.202645 0.294806 0.21554 0.283794 0.241777 0.232855 0.463198 0.234144 0.15142 0.312594 0.081076 0.220599 0.281988 0.26158 0.21031 0.258186 0.106611 0.157161 0.24907 0.300621 0.239838 0.241168 0.187459 0.33851 0.200549 0.280677 0.233996 0.295504 0.178446 0.274467 0.18351 0.272267 0.245089 0.30548 0.226976 0.319164 0.352498 0.272045 0.411657 0.317062 0.245679 0.249905 0.324815 0.212608 0.19065 0.228242 0.319123 0.161796 0.222453 0.250633 0.159693 0.188686 0.342132 0.281646 0.334669 0.089793 0.259196 0.273596 0.345855 0.289007 0.383143 0.156662 0.31944 0.234978 0.365238 0.362004 0.271136 0.31443 0.327432 0.310543 0.318961 0.327258 0.273259 0.291889 0.240509 0.323329 0.315887 0.265692 0.367584 0.265435 0.311625 0.325877 0.118957 0.352088 0.328152 0.380891 0.285197 0.272198 0.362903 0.174778 0.187733 0.367137 0.389172 0.216745 0.257631 0.208926 0.201978 0.166519 0.394435 0.203023 0.138376 0.309235 0.288946 0.232818 0.407359 0.174608 0.29806 0.319065 0.210996 0.269264 0.253151 0.270122 0.200281 0.16621 0.264865 0.24415 0.253555 0.22509 0.16207 0.226025 0.268442 0.269452 0.19971 0.35423 0.306264 0.379176 0.212138 0.360324 0.260861 0.251174 0.380434 0.173883 0.183366 0.270164 0.220242 0.224857 0.319678 0.265396 0.213789 0.257142 0.413658 0.310209 0.250293 0.185806 0.277591 0.291807 0.250875 0.221382 0.162722 0.322992 0.284653 0.360312 0.173597 0.185779 0.260096 0.217605 0.342317 0.295417 0.306089 0.28615 0.144634 0.193834 0.259994 0.306451 0.238475 0.376848 0.317684 0.317326 0.317616 0.248689 0.243321 0.267448 0.128862 0.205479 0.22382 0.443122 0.302049 0.296842 0.3343 0.273145 0.328508 0.321244 0.279993 0.239374 0.31767 0.34008 0.42074 0.248336 0.311358 0.090746 0.244562 0.221122 0.23128 0.319185 0.244818 0.318622 0.055434 0.297108 0.305041 0.233084 0.279998 0.203625 0.358332 0.373194 0.226709 0.353974 0.320693 0.344829 0.210927 0.348568 0.292217 0.27214 0.156461 0.277307 0.197162 0.169494 0.234971 0.365146 0.256547 0.314773 0.206562 0.271107 0.321583 0.300027 0.354164 0.315897 0.153527 0.23706 0.219551 0.196176 0.257843 0.19848 0.259093 0.212337 0.2907 0.204702 0.085703 0.409945 0.117321 0.372536 0.282131 0.314824 0.345557 0.198942 0.308 0.316552 0.307774 0.241196 0.186512 0.412596 0.177018 0.28145 0.261307 0.234041 0.140856 0.25519 0.390982 0.180543 0.277732 0.335294 0.193237 0.252199 0.298091 0.193329 0.249925 0.057563 0.170051 0.330266 0.379467 0.243315 0.240575 0.389301 0.192193 0.092742 0.353249 0.322236 0.247602 0.245593 0.208722 0.244553 0.276348 0.200737 0.224606 0.376989 0.201767 0.2187 0.283648 0.317926 0.217582 0.314684 0.240584 0.297079 0.246634 0.25842 0.14068 0.291048 0.277858 0.303067 0.176676 0.351598 0.232455 0.354829 0.15342 0.193008 0.25869 0.168002 0.358689 0.206944 0.283095 0.228079 0.329994 0.266831 0.241389 0.208554 0.24696 0.156228 0.241968 0.114897 0.2753 0.268019 0.316356 0.251428 0.077423 0.28525 0.378735 0.198843 0.362631 0.366286 0.209999 0.331072 0.221995 0.323125 0.275239 0.315096 0.240005 0.263246 0.200506 0.209304 0.271093 0.24257 0.106093 0.154732 0.268824 0.347673 0.298718 0.160256 0.199619 0.281827 0.252471 0.228201 0.168023 0.178034 0.239046 0.197133 0.176046 0.363498 0.238077 0.326785 0.26029 0.218455 0.226381 0.206722 0.34503 0.290004 0.253655 0.181869 0.232676 0.135437 0.328947 0.320924 0.329285 0.27907 0.075509 0.315175 0.171692 0.209655 0.203714 0.16696 0.351023 0.237497 0.407823 0.346424 0.225227 0.303389 0.300147 0.302502 0.167579 0.155249 0.22526 0.212389 0.195561 0.321324 0.173153 0.20634 0.214044 0.21151 0.252707 0.19453 0.239438 0.317107 0.1058 0.277933 0.371312 0.395137 0.347857 0.292492 0.243573 0.235935 0.220924 0.296878 0.181623 0.198708 0.309077 0.238375 0.342543 0.258484 0.320669 0.255076 0.310499 0.296306 0.320798 0.221418 0.241488 0.302548 0.289508 0.279572 0.295063 0.207934 0.452065 0.295878 0.345433 0.147496 0.335002 0.259448 0.337932 0.151333 0.360151 0.344411 0.226303 0.298183 0.306819 0.246017 0.298072 0.33209 0.238928 0.182459 0.301487 0.275405 0.211085 0.288062 0.335839 0.240441 0.193005 0.302005 0.231519 0.35293 0.231224 0.30763 0.220013 0.316056 0.251132 0.275568 0.27514 0.201922 0.123765 0.402274 0.242015 0.384441 0.271221 0.25216 0.332672 0.285663 0.372342 0.224825 0.235547 0.295728 0.108763 0.321852 0.287818 0.381859 0.251527 0.341336 0.196716 0.288388 0.221417 0.318837 0.306761 0.238302 0.227109 0.24513 0.169065 0.265219 0.217899 0.103456 0.281321 0.16332 0.274875 0.252363 0.190833 0.214859 0.339466 0.14475 0.138736 0.123124 0.253049 0.281986 0.231643 0.289638 0.307564 0.283562 0.331114 0.235401 0.292894 0.263396 0.170371 0.31145 0.354642 0.172497 0.211499 0.307332 0.30095 0.31563 0.23364 0.309406 0.302065 0.264642 0.230294 0.360702 0.288913 0.292967 0.282992 0.131652 0.308677 0.212837 0.28879 0.245339 0.343219 0.263852 0.265962 0.133104 0.285131 0.237903 0.319864 0.238229 0.150973 0.284689 0.234584 0.319902 0.164342 0.242229 0.231901 0.211714 0.286572 0.275452 0.230947 0.160158 0.287215 0.314083 0.318609 0.261333 0.272518 0.131545 0.179749 0.328615 0.170006 0.30967 0.164024 0.259255 0.16679 0.264392 0.383925 0.432548 0.260026 0.22645 0.343881 0.311708 0.255848 0.278455 0.17303 0.294584 0.224891 0.333172 0.170166 0.188996 0.183609 0.159077 0.253739 0.314755 0.115836 0.207029 0.202147 0.289098 0.34239 0.345722 0.114379 0.251842 0.179346 0.264441 0.28935 0.183967 0.246945 0.165885 0.334738 0.163771 0.206471 0.240317 0.256187 0.229578 0.203627 0.351677 0.180396 0.313365 0.280009 0.339997 0.329028 0.22445 0.270013 0.253774 0.322327 0.301919 0.222277 0.403366 0.266934 0.297115 0.295771 0.233936 0.149582 0.216102 0.324937 0.24094 0.185224 0.190376 0.196354 0.256735 0.276962 0.171561 0.261695 0.269895 0.210425 0.235471 0.289197 0.167384 0.344503 0.254137 0.241971 0.238293 0.189452 0.373158 0.209392 0.225708 0.282154 0.210567 0.276306 0.201758 0.170622 0.221923 0.239347 0.344946 0.224187 0.218236 0.250051 0.273595 0.329786 0.325256 0.257791 0.333335 0.295368 0.394522 0.300352 0.323523 0.249517 0.270865 0.198834 0.245211 0.207448 0.223647 0.35399 0.289236 0.256261 0.319033 0.227137 0.252115 0.324349 0.162225 0.223056 0.121574 0.284598 0.232363 0.330205 0.215692 0.313182 0.380527 0.259729 0.303708 0.20282 0.282258 0.236027 0.217887 0.367712 0.271558 0.267379 0.341719 0.204884 0.200059 0.335082 0.251743 0.200433 0.184918 0.240495 0.256545 0.236554 0.31041 0.215275 0.230729 0.179142 0.219456 0.280583 0.166523 0.291384 0.218031 0.195328 0.341484 0.393951 0.19367 0.212026 0.304894 0.219901 0.212505 0.173995 0.271864 0.216251 0.346598 0.228826 0.245561 0.285509 0.270702 0.282815 0.231177 0.227556 0.307301 0.155735 0.344213 0.22937 0.195784 0.135341 0.287282 0.116771 0.262321 0.324137 0.261928 0.112158 0.361812 0.188105 0.321323 0.172408 0.193537 0.283319 0.220687 0.266344 0.269757 0.249341 0.100082 0.18317 0.300605 0.304549 0.258246 0.222234 0.293187 0.326822 0.296121 0.296075 0.262126 0.335798 0.387502 0.336422 0.334213 0.239107 0.180545 0.300523 0.187278 0.25899 0.284577 0.291006 0.157797 0.293878 0.275181 0.315566 0.218473 0.252588 0.200741 0.235605 0.310268 0.323877 0.333831 0.194114 0.237324 0.332521 0.243561 0.213893 0.20986 0.139242 0.209285 0.187387 0.402265 0.413184 0.328995 0.33037 0.295258 0.258137 0.27093 0.241467 0.197348 0.245296 0.442146 0.236388 0.302205 0.187656 0.26584 0.35311 0.231767 0.134285 0.327458 0.112293 0.165703 0.349346 0.275381 0.335419 0.241331 0.308049 0.204767 0.28132 0.240202 0.095486 0.249091 0.27532 0.084164 0.29536 0.260339 0.291991 0.17085 0.34789 0.252867 0.337038 0.345668 0.288609 0.28262 0.182713 0.236742 0.325094 0.285949 0.214278 0.148075 0.224103 0.389611 0.06923 0.193032 0.365505 0.28404 0.146384 0.308144 0.197619 0.303144 0.325061 0.275841 0.268377 0.304512 0.260552 0.212902 0.248148 0.335736 0.212924 0.346276 0.354203 0.287467 0.269705 0.220741 0.299066 0.323233 0.226127 0.363112 0.341314 0.309951 0.16238 0.284349 0.279306 0.166288 0.376415 0.237524 0.301081 0.138507 0.337373 0.33524 0.22889 0.295421 0.348777 0.330663 0.326729 0.314211 0.277524 0.195472 0.252961 0.247508 0.207546 0.32565 0.360851 0.280099 0.274245 0.286094 0.195474 0.308497 0.294778 0.179554 0.240535 0.281355 0.242496 0.220731 0.291467 0.304686 0.32515 0.295867 0.34285 0.192193 0.144241 0.273342 0.198942 0.201102 0.165858 0.215757 0.309997 0.382799 0.295303 0.352471 0.204127 0.210533 0.084925 0.307589 0.183774 0.188246 0.288326 0.329967 0.127807 0.158807 0.324285 0.172129 0.296816 0.338734 0.193678 0.318178 0.30578 0.204493 0.300636 0.121109 0.277342 0.345452 0.335665 0.332633 0.357656 0.22872 0.33177 0.125575 0.226207 0.314396 0.228892 0.25618 0.308485 0.304927 0.242702 0.291184 0.307607 0.288353 0.225121 0.28612 0.212849 0.31551 0.268307 0.343238 0.291562 0.226374 0.218772 0.336845 0.232916 0.349717 0.246879 0.185218 0.219736 0.265777 0.23731 0.271197 0.247945 0.288728 0.305088 0.242123 0.205142 0.323846 0.273191 0.200856 0.233669 0.303288 0.167788 0.283257 0.286716 0.196605 0.207673 0.088664 0.215858 0.300215 0.316797 0.220993 0.103834 0.323128 0.192074 0.335513 0.198855 0.286254 0.253017 0.253015 0.283838 0.286779 0.182723 0.111313 0.225468 0.237833 0.333391 0.289807 0.302769 0.377196 0.218413 0.120398 0.159732 0.248491 0.166947 0.287561 0.215651 0.250719 0.293862 0.148887 0.387022 0.106462 0.306485 0.349328 0.216789 0.127228 0.379327 0.155982 0.286051 0.171928 0.135649 0.220047 0.25796 0.237049 0.27332 0.207019 0.172559 0.24939 0.397942 0.343593 0.298109 0.277449 0.204092 0.180399 0.267126 0.286774 0.254466 0.327236 0.29617 0.215248 0.124381 0.280948 0.2307 0.278454 0.34779 0.285127 0.30948 0.275783 0.233237 0.280896 0.238886 0.273087 0.20527 0.274891 0.146301 0.291159 0.212397 0.155014 0.320631 0.345663 0.240532 0.280662 0.311624 0.199028 0.224396 0.315272 0.164441 0.308045 0.211994 0.15123 0.31284 0.338226 0.262176 0.299695 0.299245 0.318145 0.331989 0.130962 0.191489 0.229329 0.355872 0.12754 0.284225 0.293841 0.287808 0.292728 0.222238 0.272282 0.267815 0.301005 0.282383 0.238224 0.211687 0.273789 0.314137 0.266397 0.30473 0.279014 0.290852 0.242924 0.158705 0.240364 0.164784 0.264743 0.316193 0.230892 0.318745 0.190221 0.305744 0.323368 0.347363 0.279166 0.284967 0.227707 0.35977 0.236898 0.320235 0.322143 0.249022 0.25937 0.322853 0.139159 0.349512 0.131425 0.225326 0.260353 0.318999 0.214854 0.350816 0.276384 0.297704 0.225113 0.283429 0.337736 0.271214 0.383346 0.305051 0.29651 0.22227 0.183724 0.306847 0.234977 0.356648 0.248076 0.14666 0.228214 0.266234 0.20263 0.301292 0.298122 0.296463 0.247338 0.267823 0.126387 0.26075 0.277338 0.328043 0.279623 0.160822 0.301768 0.216332 0.20847 0.244073 0.215263 0.246667 0.393496 0.378802 0.209171 0.301623 0.370794 0.196788 0.373279 0.20927 0.343028 0.237422 0.321775 0.270325 0.272652 0.202409 0.235256 0.231816 0.250599 0.043947 0.378032 0.212793 0.335259 0.390518 0.155918 0.16895 0.311677 0.344826 0.270375 0.332767 0.257774 0.385048 0.263662 0.202193 0.338831 0.299888 0.267342 0.186372 0.163847 0.348022 0.329954 0.30735 0.284122 0.179998 0.199564 0.157954 0.310039 0.218848 0.175895 0.207791 0.267007 0.219014 0.359753 0.119797 0.358697 0.262084 0.342221 0.314262 0.30069 0.249398 0.29788 0.159888 0.275953 0.287198 0.13175 0.222674 0.213362 0.262194 0.20844 0.227489 0.206092 0.244111 0.345046 0.299006 0.267145 0.306676 0.305787 0.327381 0.26461 0.310219 0.301474 0.273932 0.290544 0.162041 0.297155 0.390594 0.334139 0.21401 0.28215 0.083592 0.262953 0.29249 0.31084 0.371357 0.30625 0.181003 0.16216 0.360584 0.230359 0.238816 0.314143 0.249009 0.275362 0.261923 0.248826 0.295921 0.294404 0.278345 0.27694 0.317435 0.292897 0.37656 0.216481 0.102934 0.301957 0.30816 0.28141 0.145676 0.323394 0.208284 0.155564 0.125027 0.2805 0.375768 0.335852 0.207236 0.132681 0.222419 0.108457 0.290646 0.298033 0.255191 0.261412 0.295086 0.338001 0.318899 0.229254 0.213198 0.182047 0.321773 0.28288 0.328379 0.402485 0.340058 0.321569 0.28972 0.268193 0.238656 0.288122 0.23834 0.215611 0.225861 0.199368 0.279618 0.318925 0.252413 0.152146 0.242117 0.096843 0.321634 0.28443 0.212712 0.31651 0.186796 0.166299 0.211993 0.318955 0.308206 0.187508 0.293057 0.271422 0.231697 0.314454 0.374434 0.254917 0.137129 0.196359 0.420605 0.085946 0.169969 0.326907 0.245621 0.346324 0.324818 0.253102 0.272626 0.174725 0.272454 0.290905 0.261252 0.225346 0.346812 0.13177 0.318497 0.343358 0.27117 0.332816 0.33941 0.332522 0.141734 0.264166 0.119587 0.379674 0.288642 0.311523 0.163823 0.14533 0.181392 0.239798 0.149729 0.266502 0.178407 0.239502 0.159702 0.12362 0.279133 0.182362 0.251795 0.337515 0.239934 0.324864 0.141693 0.277499 0.290912 0.241869 0.22169 0.296949 0.235983 0.148896 0.196794 0.314612 0.305859 0.333378 0.169383 0.277378 0.344862 0.210415 0.216688 0.239554 0.17861 0.304543 0.226151 0.211202 0.257722 0.324127 0.171521 0.236602 0.305386 0.385796 0.378775 0.328948 0.257597 0.257539 0.319127 0.239276 0.093093 0.430414 0.269111 0.296466 0.172187 0.149468 0.199297 0.286833 0.186445 0.281548 0.227126 0.222579 0.182897 0.278063 0.19252 0.287834 0.306225 0.179331 0.262882 0.267687 0.385247 0.171997 0.263738 0.225598 0.259139 0.169626 0.184406 0.185561 0.23911 0.225406 0.278644 0.304367 0.146114 0.275424 0.223372 0.302702 0.198713 0.263144 0.325083 0.369972 0.254263 0.350527 0.348498 0.197817 0.363396 0.247663 0.191698 0.176476 0.171336 0.090718 0.296363 0.197119 0.248535 0.261435 0.342627 0.291793 0.257543 0.250611 0.173205 0.251781 0.337329 0.192628 0.268677 0.08561 0.341382 0.313058 0.291071 0.259597 0.321946 0.295769 0.251975 0.252868 0.237055 0.285477 0.145286 0.208402 0.32689 0.17944 0.342887 0.202395 0.239001 0.251434 0.254879 0.226907 0.208259 0.281152 0.205321 0.221197 0.215131 0.171783 0.24878 0.155366 0.202906 0.35154 0.232405 0.24663 0.256533 0.262825 0.294817 0.371971 0.318159 0.191051 0.214948 0.219913 0.233415 0.222239 0.261342 0.297664 0.392374 0.307177 0.248505 0.239041 0.292786 0.309203 0.378354 0.117169 0.244655 0.347079 0.240646 0.300053 0.130558 0.271397 0.277237 0.304445 0.187426 0.088905 0.245182 0.244571 0.362953 0.289836 0.341953 0.284485 0.340746 0.204324 0.320922 0.199503 0.286884 0.235253 0.254332 0.169032 0.223902 0.269587 0.218688 0.260147 0.228524 0.321625 0.189907 0.212029 0.378708 0.324261 0.194205 0.252551 0.417105 0.179362 0.219417 0.302233 0.27044 0.157315 0.377989 0.237709 0.19304 0.213053 0.202094 0.207985 0.2448 0.35185 0.185715 0.274096 0.283615 0.304546 0.182532 0.134724 0.283909 0.266621 0.284961 0.260339 0.26182 0.288909 0.125584 0.270113 0.243435 0.318341 0.258058 0.23857 0.207348 0.20175 0.319034 0.229303 0.27453 0.312314 0.296285 0.236041 0.257475 0.268887 0.248604 0.123434 0.149396 0.275938 0.2369 0.272057 0.25819 0.241908 0.208721 0.215853 0.237359 0.304752 0.245735 0.355451 0.321755 0.281427 0.355444 0.306951 0.338606 0.332662 0.290473 0.243263 0.211615 0.288891 0.21483 0.36244 0.286541 0.304237 0.105462 0.221477 0.14657 0.200747 0.250392 0.383929 0.255804 0.310182 0.194111 0.268262 0.271858 0.240886 0.189202 0.317554 0.232919 0.23504 0.239061 0.232781 0.202682 0.381317 0.205718 0.232964 0.240497 0.157235 0.311281 0.255084 0.238467 0.239263 0.300337 0.378016 0.330448 0.363531 0.306813 0.276604 0.177339 0.248514 0.275582 0.278847 0.364117 0.22913 0.211859 0.14287 0.366185 0.169523 0.192319 0.301056 0.218803 0.304125 0.282757 0.254124 0.2058 0.358663 0.31303 0.340634 0.131275 0.358442 0.283298 0.182804 0.251864 0.302049 0.407412 0.349828 0.259057 0.233543 0.328837 0.252346 0.272963 0.21758 0.312498 0.157082 0.371576 0.284691 0.293755 0.12981 0.270704 0.290285 0.189855 0.229113 0.093034 0.317964 0.323247 0.271535 0.203371 0.178593 0.354762 0.250112 0.325937 0.389241 0.263429 0.328759 0.218043 0.33155 0.255677 0.350977 0.305055 0.218666 0.246902 0.34151 0.217159 0.11701 0.362568 0.284045 0.372054 0.32528 0.20591 0.285278 0.243341 0.181683 0.285898 0.278885 0.281914 0.255359 0.279259 0.16145 0.364799 0.285313 0.135297 0.317774 0.335147 0.385686 0.355709 0.275845 0.261336 0.177839 0.352306 0.370152 0.200266 0.31633 0.174045 0.328202 0.106222 0.224644 0.234067 0.302758 0.383746 0.210434 0.161676 0.396917 0.153764 0.289528 0.262269 0.405803 0.309993 0.303217 0.331837 0.308763 0.251269 0.242344 0.373773 0.190001 0.298993 0.233562 0.217498 0.386616 0.282086 0.314265 0.30549 0.266 0.232916 0.368302 0.172378 0.309242 0.285309 0.274331 0.200558 0.264578 0.30082 0.117814 0.343636 0.226731 0.210872 0.369662 0.274323 0.225276 0.164034 0.256813 0.274305 0.340775 0.255318 0.310192 0.330646 0.318259 0.186365 0.178718 0.234501 0.265778 0.230132 0.369138 0.3816 0.245894 0.265204 0.140785 0.231053 0.326078 0.32268 0.280086 0.161277 0.157752 0.290387 0.233026 0.27516 0.306353 0.287267 0.235848 0.353543 0.132936 0.22438 0.301252 0.094805 0.345368 0.147356 0.155656 0.265445 0.220489 0.341035 0.226646 0.220812 0.221693 0.358683 0.265927 0.260855 0.228121 0.324378 0.248348 0.303285 0.277198 0.262312 0.202976 0.311964 0.253682 0.252598 0.268078 0.37354 0.258404 0.264703 0.193936 0.368455 0.259703 0.291748 0.262869 0.211116 0.304728 0.244516 0.355595 0.204806 0.264253 0.339203 0.340215 0.266293 0.228802 0.202811 0.356541 0.311624 0.132634 0.188688 0.381442 0.208111 0.335417 0.288985 0.18763 0.347553 0.18021 0.23189 0.19417 0.268074 0.209298 0.185087 0.257929 0.175301 0.300016 0.179602 0.314366 0.354656 0.183144 0.30245 0.284967 0.259275 0.239207 0.322867 0.33753 0.234924 0.324414 0.34763 0.364724 0.220419 0.238208 0.220716 0.271525 0.217138 0.325014 0.22057 0.250496 0.206298 0.223688 0.279608 0.310692 0.315605 0.198543 0.109146 0.292861 0.200255 0.269632 0.251764 0.163965 0.281304 0.422857 0.28019 0.37283 0.203564 0.350615 0.288467 0.275575 0.289402 0.261445 0.27333 0.214133 0.264377 0.288289 0.2683 0.262523 0.260884 0.367515 0.214209 0.254443 0.238952 0.322203 0.306603 0.264968 0.185783 0.182772 0.160846 0.314502 0.254387 0.23475 0.285461 0.17561 0.271049 0.252224 0.270212 0.222459 0.330134 0.16786 0.310916 0.200493 0.276523 0.215536 0.198609 0.125018 0.29309 0.166209 0.30457 0.315673 0.062979 0.27272 0.270356 0.28799 0.083812 0.336531 0.177504 0.22945 0.254828 0.38244 0.274905 0.308355 0.201407 0.253392 0.164006 0.274139 0.263909 0.236348 0.291119 0.2644 0.273125 0.149605 0.210374 0.227342 0.221187 0.251017 0.200999 0.221628 0.279973 0.35643 0.364198 0.16972 0.219549 0.229944 0.303827 0.177627 0.198622 0.241916 0.243828 0.179181 0.197999 0.208961 0.215487 0.139978 0.226452 0.260996 0.322363 0.252555 0.227 0.303062 0.183603 0.165293 0.230709 0.284826 0.239213 0.297567 0.397408 0.140496 0.131926 0.314908 0.304759 0.258583 0.263126 0.345949 0.296102 0.290628 0.055094 0.258746 0.209805 0.337336 0.236196 0.279348 0.297951 0.202917 0.234994 0.340196 0.069695 0.362555 0.209299 0.257622 0.14912 0.284263 0.261129 0.196164 0.328827 0.209699 0.333641 0.275364 0.212022 0.183055 0.360139 0.33945 0.194189 0.257606 0.290473 0.362025 0.218372 0.200655 0.331143 0.289294 0.325379 0.272398 0.328637 0.130736 0.296358 0.26032 0.240112 0.357234 0.246392 0.264312 0.168617 0.293322 0.11444 0.323696 0.265921 0.303242 0.244459 0.293161 0.196495 0.285008 0.328481 0.343571 0.345462 0.249663 0.333468 0.185443 0.15982 0.234316 0.293692 0.293552 0.175696 0.163987 0.366195 0.294091 0.171411 0.293774 0.246796 0.312266 0.266793 0.187651 0.293596 0.164344 0.341127 0.289235 0.289385 0.222847 0.378267 0.336288 0.305774 0.209419 0.261894 0.256519 0.365283 0.301913 0.227394 0.333416 0.23384 0.172402 0.108621 0.360855 0.216409 0.293204 0.250812 0.215286 0.180413 0.272256 0.27112 0.288615 0.222673 0.147446 0.276873 0.328091 0.252083 0.144745 0.189488 0.21932 0.236026 0.32437 0.264799 0.142223 0.346992 0.199666 0.224167 0.225058 0.150353 0.211433 0.381291 0.338767 0.346963 0.301181 0.23037 0.150378 0.218677 0.33895 0.130742 0.311762 0.125815 0.293553 0.271649 0.400887 0.228111 0.170648 0.31696 0.243847 0.256916 0.361909 0.178867 0.338376 0.297541 0.251959 0.356059 0.287227 0.281277 0.264166 0.222838 0.186517 0.307166 0.332432 0.326896 0.33027 0.299836 0.200248 0.281178 0.360207 0.225794 0.172577 0.332457 0.19634 0.290963 0.183401 0.259006 0.235802 0.240892 0.254318 0.363532 0.308782 0.23267 0.296546 0.301676 0.259621 0.153648 0.156317 0.281949 0.324668 0.30337 0.299584 0.199142 0.250579 0.101555 0.17099 0.305584 0.217499 0.241823 0.291655 0.252533 0.253157 0.190618 0.218575 0.320233 0.2178 0.234253 0.246056 0.110504 0.306595 0.271822 0.192068 0.262844 0.160987 0.342083 0.274675 0.346921 0.231241 0.219768 0.359089 0.351078 0.321532 0.255186 0.058264 0.237567 0.233327 0.296565 0.24636 0.203124 0.286678 0.370138 0.259392 0.245913 0.310188 0.314245 0.31484 0.280327 0.291422 0.225836 0.271238 0.318193 0.291681 0.264227 0.295342 0.216439 0.288615 0.240639 0.17085 0.272102 0.166145 0.326458 0.256933 0.239865 0.147665 0.212342 0.118437 0.230592 0.209192 0.270974 0.240697 0.301088 0.254485 0.231667 0.220115 0.239625 0.261075 0.323149 0.372223 0.328959 0.259911 0.247922 0.254352 0.337955 0.307044 0.272656 0.268687 0.278581 0.179439 0.197875 0.253879 0.330447 0.328451 0.322889 0.212974 0.209379 0.374132 0.342242 0.143362 0.307252 0.271674 0.206574 0.222326 0.185833 0.282052 0.274079 0.286305 0.335183 0.302398 0.260761 0.327531 0.305202 0.228644 0.266955 0.335233 0.203027 0.258739 0.151266 0.285775 0.147043 0.097069 0.270162 0.26302 0.198417 0.327864 0.200231 0.298887 0.261934 0.271381 0.130179 0.281111 0.085067 0.350975 0.295503 0.293983 0.352341 0.233062 0.262069 0.254089 0.110528 0.31781 0.229804 0.310013 0.196477 0.335331 0.228338 0.302046 0.247828 0.276125 0.179156 0.331404 0.326693 0.225013 0.22811 0.321111 0.294535 0.274964 0.291224 0.24001 0.395564 0.168361 0.241874 0.295304 0.250275 0.314293 0.095821 0.197437 0.220382 0.246412 0.265848 0.264121 0.229078 0.293026 0.192619 0.303358 0.155066 0.359441 0.218935 0.248852 0.27762 0.274458 0.279427 0.144112 0.279807 0.298004 0.144983 0.228834 0.301394 0.292062 0.141381 0.267622 0.272955 0.251278 0.361609 0.258362 0.159349 0.165878 0.315884 0.302567 0.154169 0.201111 0.224126 0.202394 0.260581 0.243367 0.148758 0.334978 0.207737 0.328283 0.234917 0.180658 0.131017 0.247021 0.174199 0.253844 0.243098 0.177866 0.25046 0.276324 0.163132 0.193785 0.259594 0.248292 0.261556 0.397782 0.216588 0.28451 0.333389 0.188347 0.278865 0.249484 0.080981 0.305684 0.165105 0.196208 0.235349 0.313135 0.296908 0.280144 0.27719 0.232022 0.317785 0.285181 0.312674 0.281279 0.153948 0.377444 0.323962 0.297582 0.351415 0.169213 0.234556 0.084647 0.322534 0.305471 0.17963 0.258954 0.286116 0.171565 0.343592 0.277123 0.227151 0.209304 0.306822 0.331132 0.208254 0.098696 0.276906 0.271294 0.297967 0.263106 0.171546 0.198274 0.3343 0.289432 0.320025 0.247293 0.261398 0.243121 0.259037 0.239776 0.336406 0.326871 0.381693 0.35712 0.263107 0.171475 0.211137 0.343201 0.292731 0.281432 0.166786 0.176065 0.342949 0.301047 0.273679 0.23997 0.203959 0.212933 0.163062 0.34481 0.204166 0.3109 0.244853 0.282646 0.244829 0.262078 0.172119 0.191095 0.24425 0.245504 0.246727 0.308953 0.269115 0.277694 0.25498 0.324467 0.265022 0.305543 0.156729 0.288568 0.29236 0.321635 0.310916 0.225164 0.279464 0.142742 0.284678 0.260484 0.28331 0.306911 0.224034 0.389355 0.129019 0.367504 0.287442 0.203251 0.193025 0.161769 0.251951 0.346063 0.320824 0.281347 0.288071 0.265628 0.208715 0.089596 0.220117 0.273168 0.249366 0.220732 0.335357 0.212171 0.303394 0.150665 0.213519 0.201592 0.147511 0.297763 0.281978 0.309211 0.315832 0.245685 0.288559 0.287948 0.315739 0.228264 0.301164 0.254752 0.357278 0.279627 0.344811 0.112783 0.205482 0.237613 0.280379 0.312748 0.164097 0.193901 0.287216 0.225945 0.292739 0.293525 0.258282 0.370381 0.181405 0.216557 0.119557 0.201131 0.325776 0.328777 0.290405 0.238527 0.36389 0.346353 0.291321 0.321291 0.173038 0.308487 0.26807 0.292092 0.190855 0.089154 0.299292 0.309151 0.251835 0.3318 0.295541 0.192742 0.264755 0.383374 0.12654 0.223426 0.277875 0.250221 0.277893 0.379053 0.199208 0.202455 0.236751 0.169953 0.293437 0.156955 0.270196 0.311168 0.303125 0.215586 0.191884 0.231015 0.327689 0.198048 0.284917 0.304499 0.31766 0.262497 0.320645 0.25295 0.234815 0.246299 0.181235 0.259995 0.255688 0.154916 0.263954 0.312366 0.280893 0.250347 0.240964 0.283456 0.257279 0.154195 0.295705 0.338138 0.284458 0.31498 0.335865 0.224132 0.277588 0.208779 0.375636 0.076258 0.199979 0.232736 0.235973 0.358709 0.27477 0.190856 0.265819 0.346553 0.274318 0.258972 0.306059 0.232667 0.355351 0.174516 0.131593 0.346057 0.199797 0.336791 0.251959 0.322071 0.255224 0.227944 0.285 0.308124 0.279135 0.349269 0.307442 0.29158 0.246251 0.281146 0.260563 0.200831 0.285185 0.378925 0.306954 0.31748 0.281383 0.249036 0.239862 0.299751 0.316148 0.189721 0.187157 0.131735 0.275221 0.29556 0.273864 0.303064 0.203343 0.301241 0.254261 0.126045 0.235276 0.26559 0.216565 0.265351 0.271086 0.134832 0.233858 0.274575 0.268286 0.261376 0.252408 0.101169 0.230876 0.355024 0.232447 0.273603 0.337041 0.29932 0.174428 0.30746 0.278793 0.220656 0.219309 0.25166 0.149244 0.221303 0.270817 0.264014 0.281065 0.096496 0.295622 0.20365 0.251258 0.307087 0.282923 0.151008 0.313166 0.280811 0.270687 0.301301 0.240688 0.278736 0.21461 0.119915 0.199606 0.262677 0.395915 0.265451 0.258914 0.213067 0.284286 0.259721 0.16812 0.263374 0.36577 0.260023 0.216278 0.286008 0.366719 0.288324 0.253298 0.248924 0.151943 0.277619 0.175362 0.217514 0.222234 0.260955 0.132892 0.255679 0.276126 0.246446 0.262599 0.195639 0.35253 0.278847 0.262468 0.264103 0.311931 0.243372 0.182137 0.215889 0.314764 0.340632 0.282769 0.310225 0.315158 0.188494 0.323527 0.235847 0.243875 0.1863 0.146463 0.217375 0.363704 0.231595 0.242009 0.29484 0.372072 0.428022 0.240435 0.283871 0.297941 0.14849 0.268803 0.298803 0.347362 0.214206 0.159522 0.241431 0.273177 0.352086 0.197902 0.204917 0.244066 0.304379 0.288954 0.337396 0.139641 0.269668 0.318648 0.346795 0.167039 0.20907 0.224457 0.361571 0.297747 0.30334 0.382257 0.262739 0.247378 0.178768 0.231244 0.226111 0.337537 0.274776 0.26992 0.268795 0.124668 0.328018 0.203489 0.24161 0.373283 0.280335 0.272706 0.292583 0.180719 0.386504 0.287526 0.182173 0.229606 0.224883 0.250481 0.160795 0.269359 0.284472 0.083252 0.332347 0.147905 0.232949 0.154568 0.22245 0.300182 0.288175 0.227341 0.18768 0.202132 0.212204 0.31064 0.27663 0.251456 0.181396 0.182031 0.156799 0.347186 0.254774 0.19814 0.324428 0.233965 0.347021 0.303012 0.357462 0.185583 0.437091 0.307696 0.15479 0.209929 0.245755 0.304691 0.26515 0.238186 0.326936 0.235965 0.29435 0.314731 0.258571 0.172967 0.261968 0.232053 0.345434 0.273055 0.212473 0.16385 0.181557 0.165971 0.225509 0.152985 0.283954 0.268233 0.345009 0.278419 0.272251 0.146445 0.284221 0.250984 0.257805 0.290486 0.249846 0.289288 0.379228 0.284907 0.158422 0.162652 0.329592 0.26945 0.077638 0.268633 0.296845 0.339312 0.299909 0.132611 0.256874 0.241339 0.289396 0.262315 0.259731 0.275869 0.230124 0.290893 0.121533 0.356797 0.309357 0.251763 0.366754 0.216542 0.188095 0.264971 0.169308 0.278168 0.293788 0.132501 0.277779 0.370249 0.344356 0.251467 0.216557 0.250259 0.208324 0.153897 0.304886 0.308336 0.269434 0.342115 0.312271 0.209312 0.320523 0.195452 0.226375 0.289549 0.259937 0.137226 0.325374 0.154258 0.178473 0.196007 0.19031 0.170832 0.269992 0.228406 0.304309 0.305953 0.181989 0.250015 0.325962 0.298923 0.222574 0.222164 0.171333 0.290648 0.225977 0.355549 0.277452 0.282252 0.229678 0.260071 0.349426 0.269344 0.235989 0.243065 0.110157 0.156521 0.254317 0.273712 0.08285 0.229298 0.210718 0.375989 0.260771 0.157839 0.291859 0.297562 0.164943 0.232546 0.235333 0.24129 0.071942 0.24531 0.298218 0.252096 0.321142 0.169445 0.231468 0.368205 0.244374 0.223818 0.262165 0.337787 0.247321 0.286867 0.371577 0.245391 0.257594 0.37125 0.152788 0.394447 0.385031 0.240557 0.263518 0.183033 0.221698 0.11521 0.350328 0.30347 0.153667 0.28602 0.343376 0.278149 0.199227 0.231007 0.288865 0.323145 0.401474 0.284552 0.219032 0.331745 0.200892 0.179822 0.306506 0.333169 0.222565 0.267777 0.273465 0.25648 0.235588 0.308104 0.315079 0.287404 0.320787 0.283204 0.23796 0.2373 0.161112 0.29989 0.310116 0.336857 0.354711 0.263076 0.272547 0.213251 0.375576 0.273492 0.270449 0.241251 0.296478 0.334562 0.334759 0.271755 0.210023 0.232249 0.163155 0.175706 0.297596 0.27115 0.342756 0.23835 0.258257 0.276199 0.203142 0.30023 0.223198 0.244319 0.321247 0.272694 0.245172 0.27111 0.343846 0.162715 0.278689 0.272418 0.298999 0.337622 0.161601 0.199994 0.249364 0.206649 0.256649 0.264993 0.277879 0.230175 0.168063 0.303693 0.12562 0.246172 0.369348 0.277891 0.158516 0.294393 0.12357 0.315751 0.179201 0.292184 0.159561 0.251155 0.29276 0.292654 0.256864 0.232206 0.359557 0.274413 0.246898 0.27027 0.230326 0.28593 0.190241 0.415802 0.272883 0.101266 0.223163 0.335146 0.31404 0.287958 0.258531 0.291392 0.282317 0.208572 0.312969 0.226823 0.200363 0.231518 0.324749 0.353641 0.165297 0.15808 0.178812 0.310757 0.192785 0.125308 0.227614 0.26811 0.185776 0.225421 0.371931 0.353575 0.206137 0.325684 0.269097 0.365637 0.325332 0.297173 0.13171 0.345968 0.22629 0.315711 0.228668 0.075711 0.122542 0.268228 0.39245 0.168704 0.317637 0.285788 0.224451 0.236705 0.265304 0.406744 0.269547 0.1544 0.088793 0.317595 0.169713 0.186162 0.267645 0.374171 0.161178 0.22111 0.187173 0.313917 0.248526 0.295571 0.350747 0.310387 0.298519 0.272359 0.258774 0.261069 0.309259 0.304006 0.253921 0.258769 0.169868 0.193863 0.256127 0.152555 0.233584 0.295023 0.362968 0.372259 0.219347 0.292185 0.344428 0.32572 0.300368 0.26062 0.275368 0.329154 0.274188 0.234662 0.349934 0.188154 0.300434 0.287307 0.256314 0.278404 0.264102 0.205336 0.224287 0.287127 0.265769 0.270967 0.311181 0.140457 0.226848 0.374234 0.11563 0.189004 0.328624 0.345673 0.35488 0.224106 0.168671 0.238604 0.277892 0.316915 0.186621 0.343728 0.316603 0.275145 0.302205 0.314505 0.285857 0.421572 0.205544 0.291614 0.255696 0.336708 0.330604 0.335889 0.291344 0.390216 0.322127 0.240428 0.326814 0.248062 0.257329 0.24518 0.179719 0.124975 0.307891 0.295129 0.308061 0.341964 0.209579 0.279643 0.137622 0.311919 0.198175 0.298281 0.256776 0.243201 0.266536 0.374998 0.159757 0.272476 0.360094 0.313812 0.173247 0.249126 0.270861 0.195526 0.208202 0.312129 0.204617 0.323202 0.204201 0.331281 0.301923 0.199674 0.294806 0.234335 0.216519 0.186156 0.264109 0.15307 0.285005 0.273883 0.320468 0.235963 0.229858 0.325715 0.272852 0.193514 0.307502 0.335956 0.272465 0.260517 0.254584 0.40772 0.232482 0.350513 0.291927 0.293961 0.310082 0.221138 0.207888 0.237902 0.293013 0.255227 0.275453 0.310035 0.347988 0.143061 0.263079 0.317412 0.328625 0.24648 0.330396 0.281672 0.231758 0.225664 0.20872 0.291828 0.302803 0.351169 0.206162 0.30973 0.269281 0.305346 0.288149 0.1796 0.333403 0.25118 0.284547 0.206019 0.2583 0.268578 0.321745 0.274716 0.216534 0.269142 0.16196 0.2789 0.156328 0.255987 0.168956 0.192363 0.27198 0.209507 0.230966 0.238787 0.164344 0.321733 0.308896 0.198359 0.241183 0.241659 0.195527 0.299745 0.297933 0.130606 0.335018 0.257245 0.254281 0.257917 0.291706 0.345436 0.227293 0.277792 0.305041 0.28776 0.211788 0.310283 0.306284 0.209466 0.287244 0.270024 0.28152 0.216599 0.242142 0.318317 0.176563 0.293596 0.226095 0.334276 0.303998 0.236405 0.322832 0.22817 0.350543 0.305423 0.336734 0.230455 0.286559 0.353819 0.364045 0.192332 0.270041 0.284966 0.225213 0.155885 0.321902 0.349836 0.302372 0.348747 0.241456 0.28292 0.220949 0.253835 0.250603 0.296786 0.142959 0.251056 0.279569 0.228584 0.122961 0.286699 0.260799 0.347241 0.200391 0.205884 0.306268 0.22944 0.266633 0.301859 0.290165 0.129258 0.354166 0.307988 0.245659 0.307015 0.308613 0.286423 0.306696 0.328331 0.282085 0.322079 0.22596 0.421013 0.096201 0.347025 0.244289 0.286519 0.173078 0.240769 0.313402 0.249768 0.326941 0.206494 0.315197 0.168315 0.343973 0.318343 0.334994 0.273407 0.253103 0.304028 0.309301 0.211632 0.166709 0.160314 0.178787 0.293648 0.183119 0.196184 0.296997 0.184313 0.273236 0.369604 0.323363 0.198963 0.210367 0.25478 0.304371 0.29245 0.278792 0.249805 0.190291 0.190966 0.304855 0.339927 0.304518 0.245348 0.29617 0.211965 0.262918 0.151617 0.368136 0.133644 0.251746 0.237225 0.241793 0.311582 0.161111 0.132604 0.185537 0.375267 0.267017 0.207912 0.291983 0.24486 0.297135 0.2867 0.184096 0.269738 0.172962 0.147864 0.332835 0.3181 0.331524 0.236316 0.412776 0.26982 0.219462 0.193813 0.292072 0.276875 0.340917 0.287848 0.25607 0.249238 0.180876 0.177888 0.175354 0.231925 0.203706 0.300431 0.260261 0.18148 0.324443 0.307168 0.291563 0.180062 0.274005 0.224367 0.373534 0.169794 0.290144 0.196143 0.241596 0.35927 0.156085 0.261523 0.306439 0.190106 0.296893 0.2903 0.308436 0.35341 0.259017 0.387425 0.294925 0.28299 0.225463 0.290431 0.240987 0.365422 0.287912 0.149667 0.269892 0.339188 0.236803 0.188763 0.259297 0.297772 0.221312 0.326449 0.290223 0.135881 0.26014 0.182085 0.142544 0.293857 0.314506 0.226559 0.256803 0.257919 0.192507 0.25023 0.181911 0.303276 0.232183 0.264846 0.283305 0.210153 0.169643 0.35531 0.102491 0.230844 0.227528 0.235697 0.390005 0.247313 0.255838 0.368466 0.15808 0.302212 0.21194 0.312451 0.307076 0.428178 0.30573 0.103877 0.337086 0.286346 0.255363 0.331837 0.228403 0.308607 0.289348 0.289817 0.209106 0.101335 0.188396 0.251 0.223746 0.3076 0.262264 0.386894 0.205336 0.212574 0.313146 0.329148 0.365098 0.209825 0.287003 0.325188 0.206085 0.241316 0.332708 0.284847 0.251173 0.29097 0.259239 0.326971 0.302261 0.336097 0.287603 0.281364 0.360547 0.348017 0.252242 0.355069 0.366245 0.291598 0.136528 0.296755 0.201612 0.260435 0.278339 0.290841 0.238216 0.233754 0.195717 0.181334 0.441457 0.312286 0.187983 0.291079 0.279104 0.220234 0.218374 0.259355 0.216462 0.227515 0.27587 0.269633 0.152073 0.140778 0.278063 0.327759 0.337933 0.380902 0.278981 0.305217 0.314195 0.272097 0.267213 0.298724 0.288468 0.233005 0.26685 0.264376 0.359226 0.306496 0.213168 0.174091 0.218606 0.281631 0.204816 0.243414 0.322111 0.282713 0.218805 0.279803 0.301693 0.267217 0.284661 0.411459 0.325064 0.285907 0.245421 0.252355 0.059922 0.284062 0.341596 0.363115 0.25409 0.287846 0.308426 0.356558 0.280342 0.258869 0.255574 0.298978 0.17086 0.155044 0.308996 0.087969 0.277143 0.335494 0.327964 0.259446 0.229912 0.250818 0.254089 0.131342 0.187413 0.197489 0.172203 0.315889 0.24844 0.353284 0.303103 0.384423 0.39793 0.218518 0.334217 0.160394 0.342552 0.322882 0.240021 0.149838 0.180642 0.243412 0.143536 0.278616 0.293869 0.267345 0.181477 0.331982 0.305192 0.279945 0.333319 0.286865 0.341938 0.199474 0.254628 0.220949 0.272414 0.293932 0.310997 0.17357 0.255555 0.302014 0.235916 0.234183 0.225875 0.283787 0.288994 0.182633 0.306425 0.279505 0.210314 0.276389 0.232066 0.135648 0.147291 0.221058 0.123483 0.223946 0.307949 0.275257 0.30462 0.213243 0.251834 0.195627 0.216576 0.3515 0.197486 0.293977 0.272005 0.33735 0.154654 0.233445 0.371901 0.181065 0.132959 0.320243 0.313943 0.314897 0.309271 0.180846 0.426909 0.146955 0.233777 0.338622 0.251838 0.251802 0.300884 0.246542 0.278253 0.271621 0.260681 0.289755 0.290461 0.27224 0.240545 0.305242 0.346765 0.272656 0.237815 0.222051 0.310204 0.320517 0.221954 0.379127 0.277986 0.371596 0.152171 0.334119 0.262842 0.360497 0.363274 0.345693 0.27168 0.370828 0.315403 0.190645 0.284238 0.151159 0.256723 0.28183 0.189672 0.380301 0.291975 0.236243 0.359856 0.275904 0.248063 0.268564 0.254618 0.26183 0.2882 0.202213 0.211488 0.298372 0.322942 0.254254 0.213238 0.169974 0.326497 0.303428 0.193333 0.359946 0.314878 0.324778 0.128973 0.317749 0.173332 0.330226 0.398597 0.132519 0.249879 0.273878 0.340881 0.317841 0.162641 0.342478 0.27434 0.159702 0.274135 0.241162 0.206335 0.214811 0.18283 0.261184 0.318443 0.220437 0.247653 0.193127 0.15192 0.151864 0.255497 0.14032 0.35305 0.265854 0.234757 0.276338 0.230901 0.233491 0.304243 0.34364 0.320279 0.274726 0.309131 0.224467 0.213159 0.25556 0.201326 0.275775 0.272794 0.247331 0.305214 0.157425 0.198181 0.313768 0.291973 0.355087 0.317368 0.232234 0.212986 0.226552 0.157187 0.295784 0.276981 0.335896 0.106441 0.344945 0.262201 0.111725 0.263474 0.303102 0.319432 0.289988 0.313177 0.243575 0.190349 0.347553 0.309739 0.260154 0.05011 0.145276 0.274721 0.267727 0.175896 0.241558 0.186776 0.285907 0.296338 0.252677 0.100521 0.231681 0.299075 0.273237 0.228486 0.326717 0.311533 0.294759 0.272374 0.146004 0.313134 0.223442 0.105285 0.30068 0.258592 0.253854 0.193279 0.204869 0.290479 0.219693 0.166744 0.202104 0.25647 0.248034 0.150172 0.19838 0.387975 0.316204 0.352022 0.177133 0.176098 0.288952 0.363378 0.241638 0.420861 0.288348 0.254872 0.285871 0.174044 0.349954 0.30681 0.387117 0.285944 0.132472 0.332511 0.343754 0.264282 0.196446 0.241292 0.24101 0.266473 0.227879 0.15225 0.261572 0.35131 0.09727 0.256646 0.293142 0.142272 0.178864 0.269846 0.197851 0.254201 0.251774 0.311277 0.278984 0.327216 0.25233 0.275498 0.174287 0.252677 0.259663 0.341269 0.226136 0.344476 0.356231 0.238225 0.266931 0.265019 0.222671 0.377902 0.244068 0.11757 0.334886 0.392504 0.217237 0.260423 0.190556 0.220758 0.272312 0.206853 0.26559 0.195131 0.261114 0.217309 0.193137 0.117476 0.312923 0.239645 0.277294 0.174762 0.30754 0.24627 0.248863 0.258348 0.22517 0.13127 0.210703 0.256152 0.312925 0.271713 0.301107 0.272857 0.218057 0.381036 0.32853 0.272532 0.182742 0.129644 0.347145 0.302972 0.317589 0.281571 0.211134 0.294433 0.357815 0.235418 0.199855 0.244317 0.182473 0.354853 0.293318 0.214159 0.163567 0.119624 0.17917 0.242896 0.265794 0.103241 0.216858 0.276551 0.273893 0.187603 0.356695 0.2122 0.235617 0.308353 0.110324 0.22677 0.166461 0.085112 0.289594 0.314947 0.207667 0.264165 0.279393 0.1636 0.208164 0.434567 0.302848 0.26863 0.258265 0.342425 0.268544 0.209905 0.277177 0.235234 0.271076 0.262214 0.183673 0.278141 0.324616 0.265052 0.136615 0.295877 0.225829 0.237909 0.218972 0.186771 0.306175 0.205836 0.115722 0.177404 0.317483 0.28169 0.335279 0.18188 0.270605 0.287702 0.20763 0.177291 0.212445 0.22767 0.268225 0.266881 0.314506 0.319072 0.305179 0.326057 0.341883 0.346091 0.284712 0.249904 0.372285 0.266141 0.157055 0.319257 0.246812 0.285981 0.154414 0.290798 0.247139 0.240378 0.234877 0.184128 0.212815 0.321918 0.24684 0.225435 0.225198 0.289285 0.213378 0.231741 0.364431 0.326728 0.312538 0.200459 0.327673 0.237034 0.320686 0.327207 0.25252 0.317701 0.228876 0.195649 0.272941 0.396658 0.218049 0.311848 0.224404 0.198599 0.225764 0.143094 0.199751 0.254308 0.230984 0.290827 0.34052 0.172731 0.246621 0.186237 0.370809 0.342998 0.323528 0.275665 0.172087 0.210931 0.348954 0.251689 0.285866 0.263159 0.365297 0.238305 0.231377 0.302055 0.329563 0.341516 0.237167 0.182557 0.33971 0.386828 0.200648 0.162873 0.30344 0.393778 0.165803 0.376166 0.223024 0.17509 0.287034 0.312333 0.273127 0.310493 0.152893 0.255664 0.387404 0.270092 0.257961 0.262816 0.259146 0.203503 0.257471 0.272375 0.337012 0.362381 0.366576 0.311333 0.222637 0.19587 0.249909 0.390765 0.221602 0.090277 0.319204 0.288551 0.278623 0.248876 0.210632 0.28585 0.246482 0.373097 0.347332 0.349486 0.219391 0.241672 0.224148 0.219121 0.338892 0.323303 0.186116 0.332969 0.357186 0.25921 0.27758 0.212478 0.267892 0.228304 0.184897 0.281762 0.299989 0.205664 0.362541 0.173758 0.346192 0.241408 0.294169 0.288394 0.29512 0.26091 0.200037 0.271704 0.290215 0.15567 0.334764 0.224245 0.330224 0.291062 0.259754 0.26158 0.28466 0.310825 0.276292 0.305987 0.254244 0.334675 0.217717 0.297947 0.312615 0.216551 0.230884 0.228348 0.261034 0.295385 0.16825 0.255663 0.292741 0.206384 0.283461 0.233313 0.19913 0.236918 0.351181 0.284751 0.314943 0.327545 0.171396 0.227985 0.285054 0.162164 0.211789 0.333912 0.186654 0.306287 0.28145 0.367233 0.223641 0.235324 0.156044 0.265379 0.22109 0.186911 0.295107 0.283629 0.209336 0.290504 0.240789 0.176943 0.24212 0.293748 0.29277 0.282093 0.369104 0.216473 0.233348 0.300196 0.239461 0.329875 0.274266 0.262126 0.387032 0.244895 0.285857 0.232231 0.258438 0.279279 0.250996 0.234139 0.128069 0.281133 0.272556 0.258586 0.197796 0.305739 0.221756 0.281319 0.259736 0.239111 0.266907 0.2598 0.192855 0.323678 0.280738 0.2366 0.266337 0.244887 0.215205 0.205689 0.221061 0.328341 0.283356 0.199859 0.283209 0.24366 0.304922 0.272045 0.277517 0.241694 0.179613 0.212174 0.217332 0.249899 0.160418 0.276501 0.214144 0.345918 0.178591 0.28192 0.2336 0.205871 0.336193 0.323692 0.233669 0.226267 0.245076 0.197645 0.069489 0.22958 0.292164 0.229401 0.289056 0.24724 0.239317 0.311384 0.299027 0.184793 0.333932 0.2531 0.267343 0.299221 0.236435 0.344383 0.231653 0.28936 0.111261 0.18803 0.239939 0.226941 0.254514 0.155648 0.222421 0.121869 0.369937 0.405906 0.42016 0.297284 0.301711 0.397903 0.337936 0.266894 0.152385 0.314234 0.147577 0.199941 0.311748 0.242502 0.29575 0.236435 0.248369 0.33572 0.201165 0.350935 0.233915 0.1531 0.304992 0.144637 0.301154 0.33264 0.359882 0.301905 0.255803 0.259737 0.342704 0.30864 0.200532 0.276478 0.357993 0.28713 0.242784 0.345288 0.331812 0.232206 0.277166 0.250976 0.305239 0.243122 0.212344 0.290804 0.13211 0.158625 0.238846 0.186129 0.265397 0.18803 0.321882 0.277986 0.250366 0.146042 0.172857 0.281872 0.246732 0.14965 0.281082 0.263017 0.325244 0.189613 0.253082 0.333005 0.341776 0.206513 0.254923 0.336767 0.301194 0.246897 0.333128 0.223807 0.348468 0.175655 0.293727 0.221405 0.223499 0.292726 0.25352 0.281478 0.103461 0.259713 0.22279 0.239024 0.260725 0.197818 0.18136 0.361663 0.242709 0.188776 0.281054 0.119874 0.274553 0.242364 0.262012 0.194359 0.329861 0.327225 0.294472 0.303176 0.280576 0.301576 0.195058 0.243802 0.226419 0.21456 0.3351 0.134459 0.366138 0.322205 0.203257 0.255079 0.257262 0.213508 0.355557 0.251396 0.238994 0.052872 0.313863 0.197003 0.200303 0.362996 0.260969 0.299402 0.28183 0.34787 0.298528 0.336034 0.122567 0.325558 0.324798 0.34086 0.234875 0.207384 0.214585 0.30574 0.17048 0.379252 0.253712 0.356424 0.196633 0.173744 0.222762 0.356334 0.279651 0.200152 0.343727 0.237662 0.319663 0.253774 0.176367 0.331888 0.235042 0.384102 0.175682 0.155091 0.176414 0.322405 0.192046 0.203942 0.250359 0.371159 0.34637 0.145271 0.282451 0.300819 0.293446 0.308666 0.29303 0.288046 0.277544 0.278354 0.292622 0.350035 0.213482 0.277548 0.292161 0.234881 0.179361 0.197762 0.3033 0.073823 0.174625 0.295316 0.085679 0.219097 0.30126 0.242851 0.218146 0.128492 0.313599 0.33777 0.274551 0.304427 0.110811 0.210232 0.24912 0.299352 0.200304 0.339141 0.181614 0.256775 0.259807 0.342806 0.133949 0.293727 0.304652 0.267489 0.272028 0.202257 0.219406 0.271417 0.313472 0.326389 0.268981 0.316748 0.224207 0.320759 0.24764 0.302151 0.273079 0.283201 0.275478 0.313075 0.313126 0.260015 0.221723 0.280377 0.296897 0.26863 0.258568 0.261709 0.222989 0.320237 0.319709 0.115018 0.297167 0.166063 0.176172 0.170016 0.173541 0.281897 0.378804 0.316634 0.198067 0.220405 0.13155 0.262762 0.377617 0.39262 0.288531 0.258788 0.243601 0.240739 0.331285 0.164671 0.323605 0.287907 0.299221 0.220317 0.317162 0.303695 0.354527 0.271654 0.319873 0.171179 0.267046 0.267899 0.107743 0.268402 0.189315 0.262295 0.154932 0.2958 0.11907 0.269906 0.168995 0.185037 0.246099 0.155391 0.285035 0.299752 0.288547 0.204407 0.154542 0.290784 0.119925 0.251439 0.323668 0.426331 0.241224 0.078749 0.351264 0.321144 0.227024 0.290322 0.261168 0.141754 0.305907 0.242657 0.302167 0.269817 0.298451 0.233034 0.312472 0.144677 0.211009 0.258183 0.254031 0.294903 0.197146 0.220414 0.291446 0.276214 0.352874 0.217822 0.233732 0.265788 0.264523 0.165381 0.280181 0.266448 0.263119 0.365993 0.114552 0.312365 0.312994 0.294177 0.342717 0.319958 0.152283 0.170999 0.183018 0.326684 0.296433 0.140561 0.278916 0.151273 0.35143 0.152729 0.307247 0.36541 0.346819 0.311928 0.174115 0.267778 0.325265 0.312938 0.354467 0.383525 0.265051 0.277398 0.202071 0.175006 0.321393 0.262277 0.225918 0.271835 0.16647 0.251654 0.287318 0.278826 0.301331 0.097201 0.263705 0.305714 0.219848 0.204015 0.246121 0.266511 0.257501 0.27057 0.126785 0.181721 0.135883 0.118716 0.219606 0.319048 0.280441 0.271174 0.410081 0.23837 0.29245 0.165957 0.198761 0.17425 0.204421 0.363854 0.278035 0.352719 0.220695 0.227228 0.220112 0.286599 0.202268 0.275529 0.296664 0.174669 0.303175 0.261295 0.285407 0.330919 0.386476 0.267782 0.125556 0.21296 0.119394 0.223276 0.326389 0.25141 0.258478 0.288896 0.250387 0.386924 0.350172 0.393162 0.19779 0.264948 0.268672 0.297558 0.316509 0.288319 0.270183 0.299884 0.356931 0.277699 0.306373 0.164877 0.15939 0.226616 0.35641 0.22774 0.232414 0.404923 0.232416 0.270283 0.392762 0.180129 0.245303 0.324437 0.245258 0.281634 0.291102 0.340611 0.299566 0.303049 0.274702 0.293891 0.283279 0.2309 0.309547 0.374801 0.326944 0.239365 0.362957 0.124427 0.31164 0.198526 0.320441 0.295916 0.358944 0.3471 0.162031 0.214955 0.192356 0.229602 0.275288 0.233065 0.385175 0.251668 0.300014 0.331204 0.24635 0.222028 0.357005 0.273582 0.162131 0.300752 0.337115 0.22669 0.323784 0.291934 0.359613 0.285018 0.202505 0.149342 0.172857 0.193196 0.203059 0.280858 0.303558 0.339656 0.240039 0.199562 0.289488 0.320551 0.224761 0.352093 0.311509 0.382887 0.22037 0.229869 0.313252 0.193309 0.197887 0.365233 0.337125 0.182951 0.219554 0.233322 0.286237 0.218428 0.304182 0.160164 0.19938 0.34608 0.209575 0.272009 0.260776 0.208589 0.120907 0.242788 0.169116 0.238133 0.283455 0.391781 0.299195 0.134656 0.21092 0.169701 0.310558 0.244849 0.090892 0.294217 0.108864 0.217132 0.357736 0.112443 0.155809 0.247567 0.226058 0.290349 0.28138 0.157451 0.14139 0.404164 0.314065 0.25017 0.270421 0.198873 0.265369 0.281678 0.368816 0.071814 0.078313 0.203903 0.236243 0.324192 0.243256 0.248103 0.195388 0.280391 0.241551 0.216743 0.235782 0.258666 0.324648 0.23552 0.417312 0.308468 0.332662 0.27885 0.117723 0.294037 0.283644 0.228379 0.335466 0.132719 0.112652 0.199596 0.352856 0.209867 0.216483 0.326819 0.267593 0.24794 0.199138 0.234434 0.402923 0.3116 0.32504 0.230959 0.348359 0.343532 0.121039 0.221973 0.284365 0.304249 0.205658 0.343548 0.311794 0.122131 0.292118 0.271954 0.328583 0.285507 0.380717 0.243557 0.264463 0.348688 0.333757 0.265786 0.310302 0.238415 0.273657 0.384317 0.275162 0.213304 0.321413 0.21008 0.228372 0.179409 0.273547 0.342813 0.275973 0.273823 0.143247 0.175441 0.252459 0.145539 0.180397 0.240591 0.312618 0.321149 0.209884 0.192851 0.214135 0.389752 0.249036 0.18823 0.285322 0.321998 0.145671 0.256089 0.274848 0.260869 0.415606 0.343626 0.198115 0.244428 0.267592 0.238798 0.229447 0.223437 0.188043 0.249679 0.249403 0.309932 0.230765 0.340508 0.238203 0.260096 0.2603 0.282449 0.260305 0.2003 0.294036 0.297338 0.349683 0.278509 0.247005 0.263433 0.233812 0.167065 0.329752 0.200135 0.114633 0.190981 0.320253 0.196068 0.217836 0.23896 0.205656 0.03953 0.287705 0.13887 0.419035 0.296984 0.258137 0.298828 0.21589 0.240899 0.221867 0.159933 0.358932 0.288535 0.269827 0.2436 0.106792 0.275128 0.264749 0.340498 0.232633 0.311849 0.19853 0.305746 0.312512 0.284949 0.242987 0.291992 0.096132 0.311395 0.17378 0.080942 0.283658 0.364953 0.221153 0.295254 0.290021 0.294323 0.290406 0.338776 0.300127 0.239546 0.174396 0.321623 0.299531 0.307359 0.275219 0.141007 0.17375 0.084827 0.198798 0.198951 0.319177 0.250854 0.19551 0.236248 0.174864 0.064112 0.30688 0.294085 0.284289 0.296828 0.212892 0.232184 0.308286 0.357291 0.297294 0.34872 0.238178 0.19848 0.148494 0.265843 0.215004 0.266195 0.300916 0.223389 0.240251 0.141673 0.2937 0.274252 0.250098 0.329972 0.211795 0.2968 0.184601 0.295581 0.312059 0.180175 0.285949 0.173683 0.182347 0.328562 0.266465 0.145945 0.341297 0.306941 0.218971 0.289316 0.076728 0.361698 0.31397 0.184358 0.295842 0.338023 0.343666 0.331279 0.27945 0.290804 0.382954 0.266182 0.166929 0.077353 0.190641 0.274353 0.340213 0.181503 0.335892 0.069281 0.25612 0.183189 0.28213 0.150581 0.233132 0.224146 0.154999 0.210032 0.201262 0.288931 0.229721 0.251862 0.193276 0.338369 0.312927 0.264284 0.298626 0.219828 0.223423 0.408165 0.229004 0.224519 0.150076 0.26663 0.131063 0.272573 0.323531 0.184944 0.217374 0.321591 0.270611 0.245582 0.363498 0.145132 0.291021 0.357298 0.258077 0.256185 0.107843 0.358389 0.310195 0.248376 0.380221 0.24444 0.236362 0.241013 0.345911 0.297978 0.089979 0.261145 0.204423 0.238572 0.247371 0.285238 0.354853 0.285803 0.33123 0.32063 0.190188 0.227571 0.311597 0.280869 0.347773 0.256726 0.236351 0.314039 0.36318 0.264681 0.226133 0.297156 0.128342 0.201717 0.077443 0.242877 0.177224 0.222239 0.262224 0.122218 0.249704 0.260567 0.236542 0.218869 0.295633 0.27035 0.2797 0.211238 0.201719 0.229908 0.349011 0.20316 0.224411 0.279821 0.30423 0.222527 0.216602 0.261755 0.289767 0.24857 0.226359 0.278455 0.289436 0.154015 0.261749 0.291496 0.200541 0.26505 0.281483 0.279441 0.180785 0.210248 0.408917 0.276216 0.207754 0.218931 0.301895 0.366769 0.208411 0.332902 0.2403 0.340195 0.179177 0.323321 0.339571 0.290567 0.289843 0.298041 0.33677 0.33134 0.232202 0.208136 0.35008 0.183848 0.142069 0.341435 0.231347 0.316668 0.261526 0.340857 0.209644 0.287785 0.247579 0.171975 0.199841 0.298407 0.249349 0.263967 0.219018 0.270932 0.299275 0.343197 0.221044 0.270001 0.257344 0.16623 0.335278 0.151651 0.197209 0.35793 0.31316 0.196441 0.084084 0.306939 0.283949 0.241612 0.172253 0.254857 0.198386 0.326249 0.277065 0.283114 0.2663 0.097086 0.275999 0.260431 0.084361 0.323256 0.263493 0.263008 0.329088 0.179492 0.175431 0.346886 0.170872 0.282097 0.30829 0.225231 0.101887 0.33381 0.168691 0.282655 0.173036 0.265173 0.297987 0.150909 0.287084 0.287597 0.346969 0.193023 0.285692 0.268423 0.332244 0.254835 0.232168 0.231037 0.347463 0.205189 0.256128 0.22555 0.265471 0.294533 0.284563 0.208201 0.277811 0.266904 0.198152 0.265172 0.317274 0.322809 0.236018 0.205774 0.305579 0.303052 0.254601 0.318776 0.081494 0.258764 0.365649 0.339353 0.250929 0.197957 0.22683 0.295144 0.344886 0.229943 0.243378 0.247293 0.218122 0.287727 0.224184 0.338777 0.318506 0.322805 0.302491 0.143798 0.288192 0.209479 0.288436 0.193025 0.290026 0.19635 0.344235 0.258734 0.310063 0.296294 0.302789 0.331584 0.34249 0.216881 0.240267 0.237089 0.258628 0.26838 0.126792 0.177875 0.299852 0.337822 0.29668 0.068028 0.310954 0.285473 0.214142 0.327902 0.287497 0.361928 0.297201 0.268624 0.321162 0.289023 0.275754 0.292945 0.083613 0.321303 0.249903 0.303206 0.107031 0.23636 0.306561 0.247024 0.250455 0.325762 0.272966 0.261911 0.227125 0.20994 0.294775 0.242356 0.188333 0.347984 0.238883 0.260578 0.356842 0.092752 0.3175 0.221869 0.217789 0.311185 0.158489 0.370314 0.317075 0.115674 0.218984 0.229663 0.345789 0.241118 0.178787 0.249382 0.284837 0.196127 0.26086 0.319772 0.261611 0.253499 0.326721 0.240328 0.228431 0.34646 0.252838 0.27905 0.221356 0.239169 0.270084 0.257903 0.274394 0.239956 0.368083 0.30765 0.345226 0.320743 0.205681 0.149569 0.299849 0.303047 0.17538 0.217862 0.313669 0.275542 0.132039 0.335467 0.258749 0.256107 0.279095 0.210432 0.177704 0.2272 0.353197 0.27531 0.124843 0.259665 0.241711 0.215032 0.35756 0.27164 0.272569 0.236151 0.215933 0.320398 0.153872 0.163681 0.350844 0.194476 0.330842 0.284979 0.17059 0.323587 0.378732 0.263378 0.141332 0.251376 0.30726 0.313794 0.348476 0.167966 0.232503 0.304817 0.203386 0.172904 0.293581 0.181874 0.229688 0.272714 0.262345 0.342251 0.275204 0.203103 0.308055 0.25056 0.230025 0.174226 0.282003 0.294885 0.329303 0.322943 0.231607 0.178649 0.224808 0.345799 0.335552 0.267332 0.232756 0.331982 0.267222 0.294055 0.170101 0.249134 0.219501 0.267054 0.190449 0.383011 0.138874 0.378369 0.187571 0.17149 0.262954 0.193477 0.306998 0.374613 0.310152 0.373561 0.368253 0.362033 0.258526 0.359901 0.363516 0.314036 0.30155 0.188532 0.235122 0.220809 0.261856 0.319742 0.287425 0.159027 0.333405 0.291111 0.209603 0.263241 0.311257 0.084194 0.253174 0.201452 0.36691 0.280418 0.359128 0.241107 0.17586 0.289898 0.249306 0.309665 0.203299 0.151423 0.250125 0.159925 0.255075 0.356208 0.264119 0.159512 0.227409 0.136494 0.286455 0.251735 0.325095 0.215789 0.281856 0.251782 0.235573 0.211341 0.284299 0.25339 0.248184 0.246778 0.130036 0.303897 0.384234 0.165938 0.251462 0.292017 0.188094 0.305873 0.315894 0.273336 0.193238 0.270726 0.289994 0.237859 0.227219 0.260413 0.298166 0.245014 0.34646 0.225451 0.380727 0.304205 0.11007 0.244072 0.193908 0.190741 0.266879 0.254161 0.234747 0.304041 0.269618 0.223387 0.355627 0.177961 0.271186 0.199168 0.19774 0.271938 0.23163 0.252996 0.181495 0.206845 0.064433 0.329253 0.371507 0.326473 0.139547 0.219283 0.112726 0.239079 0.250282 0.365479 0.31859 0.18781 0.368066 0.306146 0.34539 0.259312 0.262337 0.303547 0.242599 0.267411 0.119885 0.276645 0.322712 0.234998 0.386896 0.243934 0.376508 0.243805 0.352446 0.331472 0.232551 0.240203 0.114423 0.140641 0.256656 0.219667 0.155176 0.242114 0.275486 0.226976 0.198735 0.219309 0.186996 0.140119 0.336674 0.166413 0.258266 0.342251 0.222745 0.239003 0.197831 0.312424 0.291943 0.33303 0.198658 0.350422 0.236151 0.318116 0.258775 0.190078 0.345718 0.298333 0.217385 0.344566 0.129915 0.235738 0.383763 0.245872 0.117492 0.28736 0.20313 0.280129 0.256132 0.312733 0.338438 0.33846 0.309237 0.361732 0.441653 0.235285 0.275004 0.334287 0.155436 0.327748 0.233499 0.230496 0.285828 0.318539 0.244124 0.339877 0.300116 0.088678 0.297414 0.214135 0.345012 0.120185 0.148017 0.322858 0.292557 0.093309 0.214551 0.322217 0.330066 0.306487 0.315637 0.251002 0.247673 0.246985 0.252662 0.167132 0.257809 0.354822 0.258678 0.140758 0.297708 0.08582 0.097419 0.216293 0.216078 0.289657 0.259099 0.27448 0.208862 0.224119 0.221833 0.202592 0.283052 0.245877 0.269201 0.29024 0.378493 0.29402 0.316568 0.209405 0.277946 0.374661 0.266082 0.232937 0.244737 0.245008 0.277034 0.287612 0.174773 0.159909 0.227426 0.297545 0.248256 0.260941 0.147872 0.305325 0.134437 0.21911 0.205284 0.262035 0.288465 0.289416 0.293296 0.217325 0.234312 0.377925 0.278744 0.264337 0.209602 0.316514 0.227156 0.129971 0.363526 0.237046 0.279001 0.29807 0.217313 0.104185 0.168007 0.186351 0.2041 0.206141 0.110154 0.170929 0.307527 0.157778 0.288522 0.223757 0.14347 0.283363 0.365456 0.268182 0.176781 0.288594 0.275887 0.359007 0.34473 0.281645 0.220371 0.104772 0.246395 0.236839 0.296299 0.147957 0.258101 0.202925 0.362915 0.375691 0.097244 0.426091 0.346923 0.321603 0.29314 0.287769 0.245405 0.224894 0.261829 0.178308 0.278308 0.29613 0.245696 0.24421 0.327227 0.183194 0.356097 0.410754 0.253834 0.27468 0.289292 0.221336 0.225934 0.205156 0.219185 0.295352 0.357765 0.209024 0.287691 0.250214 0.31359 0.165471 0.320909 0.366754 0.204113 0.389895 0.183655 0.238348 0.250503 0.276337 0.127921 0.250136 0.27633 0.3177 0.33999 0.232336 0.284524 0.278849 0.186564 0.383736 0.22474 0.201946 0.116477 0.197259 0.244811 0.255222 0.359394 0.402608 0.314433 0.270116 0.257539 0.316975 0.295451 0.24821 0.221677 0.215649 0.269221 0.148292 0.260801 0.351471 0.155447 0.076021 0.21073 0.24322 0.202383 0.336026 0.294074 0.249684 0.275249 0.210229 0.386192 0.20998 0.263154 0.28322 0.268414 0.141777 0.224907 0.27927 0.270032 0.294421 0.277203 0.202608 0.273729 0.426582 0.133108 0.299579 0.292112 0.179101 0.349278 0.224663 0.232429 0.107432 0.19994 0.232416 0.184332 0.309718 0.298116 0.290239 0.313661 0.151959 0.219314 0.231805 0.264665 0.281394 0.305858 0.315988 0.219598 0.25114 0.272557 0.272211 0.27507 0.26405 0.295229 0.294076 0.097313 0.267773 0.149589 0.290969 0.18278 0.335572 0.426011 0.208317 0.203324 0.214667 0.324895 0.170326 0.182972 0.30764 0.274931 0.13559 0.300573 0.221335 0.175213 0.384004 0.245656 0.277862 0.289083 0.170856 0.327901 0.338087 0.269876 0.297337 0.29382 0.224656 0.312613 0.30005 0.270179 0.229833 0.32615 0.282729 0.253318 0.210255 0.43677 0.327108 0.213436 0.174429 0.270157 0.349732 0.339721 0.271212 0.185078 0.200417 0.238688 0.241914 0.175217 0.208986 0.25054 0.276329 0.331086 0.30408 0.268809 0.280911 0.271787 0.234438 0.224444 0.261977 0.308435 0.237604 0.218901 0.303717 0.37099 0.244334 0.201124 0.324349 0.271867 0.186611 0.249573 0.318997 0.265895 0.236571 0.265088 0.274394 0.278534 0.258305 0.394682 0.209721 0.296465 0.260268 0.236395 0.364033 0.28026 0.128291 0.254038 0.228562 0.193665 0.170752 0.372404 0.372575 0.221787 0.252527 0.20949 0.328185 0.244578 0.330945 0.331424 0.292559 0.295869 0.224098 0.182035 0.318698 0.124125 0.256769 0.188336 0.316642 0.301155 0.281266 0.18561 0.2538 0.224757 0.273507 0.333744 0.236206 0.294727 0.291431 0.316449 0.219943 0.253288 0.359615 0.262515 0.222026 0.301726 0.261375 0.189552 0.234432 0.27527 0.283182 0.349989 0.206825 0.367199 0.381589 0.225922 0.153768 0.143956 0.245283 0.198281 0.242841 0.204585 0.101082 0.200362 0.176801 0.193645 0.153833 0.261368 0.289718 0.263339 0.326808 0.29932 0.246389 0.203323 0.303402 0.270291 0.231189 0.366886 0.340995 0.275213 0.271038 0.157512 0.236804 0.298975 0.199877 0.37019 0.320317 0.320004 0.31465 0.244215 0.195112 0.284285 0.291405 0.205786 0.323507 0.208496 0.320451 0.236376 0.30262 0.23102 0.28461 0.291806 0.293167 0.21684 0.147162 0.256091 0.288604 0.28748 0.23464 0.276042 0.33479 0.234539 0.303893 0.283889 0.278098 0.249293 0.348433 0.259785 0.355309 0.354792 0.206524 0.261397 0.249539 0.329575 0.35043 0.2631 0.132338 0.338512 0.325241 0.132531 0.272817 0.236422 0.357307 0.263541 0.159998 0.326031 0.179306 0.167209 0.163294 0.336213 0.126993 0.181683 0.302692 0.236393 0.360508 0.23966 0.160321 0.270659 0.27436 0.370983 0.282463 0.287589 0.275022 0.244168 0.29318 0.272254 0.196411 0.30879 0.293313 0.124659 0.332376 0.1982 0.379639 0.241304 0.204366 0.278307 0.251371 0.327881 0.303057 0.209708 0.184934 0.319846 0.307431 0.313803 0.278767 0.390985 0.305077 0.333041 0.136029 0.264536 0.26344 0.232201 0.216569 0.168497 0.258731 0.174586 0.246098 0.263294 0.221704 0.256201 0.280709 0.244479 0.24088 0.336662 0.114846 0.286423 0.243378 0.283549 0.339893 0.309062 0.339335 0.201095 0.157741 0.293209 0.296929 0.267926 0.297158 0.30697 0.425055 0.389935 0.143201 0.320601 0.232203 0.32638 0.201891 0.174471 0.309952 0.274126 0.409148 0.285599 0.211674 0.260302 0.279018 0.183855 0.379537 0.27209 0.238714 0.09381 0.223765 0.308199 0.273731 0.262086 0.27163 0.326896 0.273442 0.294775 0.321116 0.324234 0.190452 0.206518 0.237842 0.307678 0.297497 0.169533 0.295497 0.23043 0.269158 0.239453 0.305159 0.336104 0.169008 0.306227 0.340339 0.276646 0.327492 0.215231 0.291105 0.157314 0.231105 0.273605 0.198846 0.309291 0.320108 0.30594 0.206701 0.337969 0.310992 0.258735 0.248432 0.209426 0.326768 0.292325 0.31771 0.224795 0.282066 0.362464 0.300781 0.349412 0.269842 0.359972 0.280195 0.363548 0.238192 0.262291 0.287856 0.293429 0.19787 0.168585 0.214581 0.254093 0.168738 0.275017 0.180049 0.250639 0.318461 0.225384 0.053872 0.31435 0.243241 0.141393 0.228151 0.266402 0.375865 0.202173 0.203512 0.249863 0.252228 0.346494 0.283157 0.336288 0.235458 0.26996 0.174308 0.301866 0.295871 0.191386 0.399045 0.303285 0.239672 0.192984 0.231612 0.208291 0.274343 0.282892 0.242467 0.303382 0.26649 0.219477 0.327755 0.269062 0.22788 0.266465 0.347994 0.264786 0.288912 0.285953 0.247751 0.272483 0.246435 0.205876 0.396601 0.275548 0.353136 0.317644 0.179342 0.325462 0.300137 0.21613 0.263552 0.309696 0.37261 0.308051 0.235564 0.260841 0.134608 0.151603 0.345 0.252435 0.294119 0.34081 0.329106 0.232949 0.244218 0.249192 0.139702 0.269488 0.268546 0.225251 0.230773 0.267334 0.276731 0.292089 0.252294 0.32161 0.218537 0.316988 0.246466 0.316129 0.206289 0.296403 0.138075 0.181753 0.269786 0.178193 0.29619 0.380761 0.232707 0.290748 0.107187 0.144527 0.25691 0.338802 0.260036 0.275449 0.228705 0.266036 0.289198 0.29919 0.184119 0.331023 0.260627 0.166671 0.227103 0.254457 0.311868 0.183186 0.186547 0.335651 0.227566 0.279146 0.212547 0.225181 0.478447 0.188914 0.105981 0.308488 0.223175 0.340196 0.3135 0.380023 0.308805 0.322755 0.355549 0.077163 0.181931 0.218677 0.095756 0.342988 0.326355 0.242394 0.198155 0.160207 0.343103 0.279531 0.14715 0.208584 0.251667 0.171094 0.267476 0.314555 0.281024 0.333619 0.322306 0.34655 0.196399 0.25708 0.153088 0.298571 0.210395 0.274107 0.370103 0.328478 0.220696 0.269099 0.105921 0.342537 0.256895 0.264977 0.243383 0.196083 0.330926 0.392934 0.236277 0.361068 0.230535 0.203296 0.187025 0.220734 0.286652 0.305544 0.241955 0.30717 0.261209 0.080256 0.359541 0.285134 0.252101 0.297827 0.318755 0.378899 0.315528 0.318906 0.123439 0.363945 0.237216 0.28586 0.255634 0.234987 0.235676 0.271058 0.292758 0.208535 0.339947 0.187535 0.29976 0.282175 0.242705 0.312686 0.23487 0.262355 0.233345 0.335216 0.329647 0.288167 0.35385 0.152685 0.212513 0.272083 0.277736 0.180727 0.260864 0.358695 0.33972 0.190209 0.282683 0.183658 0.286803 0.304799 0.189687 0.277306 0.238838 0.271881 0.199545 0.287132 0.29746 0.303503 0.20259 0.225376 0.2389 0.212955 0.142597 0.217545 0.169001 0.284822 0.225853 0.354213 0.209634 0.398429 0.329603 0.273914 0.221091 0.220929 0.277851 0.28922 0.296693 0.281217 0.303618 0.20794 0.311314 0.342162 0.056414 0.26724 0.278638 0.150388 0.382542 0.265496 0.123452 0.297265 0.23969 0.266957 0.152126 0.219384 0.344864 0.194214 0.223215 0.150614 0.172769 0.249079 0.090471 0.210498 0.329915 0.31964 0.235452 0.323031 0.351659 0.292906 0.348321 0.318556 0.232728 0.348631 0.197153 0.189039 0.295933 0.297071 0.321919 0.310536 0.355194 0.152043 0.340186 0.146737 0.249293 0.247781 0.231817 0.337161 0.214707 0.328853 0.243717 0.379742 0.30082 0.210939 0.213101 0.282171 0.325154 0.244359 0.339726 0.200824 0.201443 0.134417 0.273068 0.279655 0.336515 0.223083 0.261554 0.223215 0.317425 0.192935 0.243238 0.169447 0.336497 0.228367 0.212858 0.329787 0.231517 0.347979 0.289589 0.229403 0.201735 0.296019 0.306802 0.325744 0.306046 0.284842 0.287992 0.286647 0.258622 0.263744 0.285205 0.23056 0.223481 0.368379 0.333681 0.189645 0.267597 0.246449 0.248407 0.348502 0.289837 0.272434 0.211389 0.202999 0.284338 0.26364 0.258201 0.299384 0.20809 0.285097 0.220675 0.191771 0.181906 0.348373 0.329263 0.241247 0.258036 0.32991 0.278163 0.150802 0.155226 0.325187 0.210276 0.219947 0.323775 0.235736 0.323557 0.297879 0.2149 0.164673 0.180234 0.223913 0.275317 0.321727 0.307287 0.210369 0.314947 0.197974 0.292526 0.27087 0.212321 0.360483 0.258732 0.299153 0.250937 0.264993 0.216913 0.263868 0.339678 0.242434 0.39992 0.285242 0.302435 0.20161 0.31963 0.2213 0.212469 0.417333 0.237276 0.193173 0.190591 0.174408 0.257224 0.257411 0.32905 0.232595 0.282296 0.105893 0.293973 0.342678 0.337301 0.232084 0.299612 0.279341 0.140128 0.247736 0.212253 0.29446 0.33481 0.247166 0.331025 0.388484 0.288964 0.239219 0.382272 0.241403 0.241898 0.180939 0.306964 0.342917 0.31934 0.276467 0.210192 0.282135 0.118955 0.182122 0.236089 0.292645 0.35453 0.252798 0.254005 0.257236 0.267142 0.251218 0.310739 0.368101 0.16673 0.335606 0.327439 0.31788 0.292753 0.264529 0.349209 0.444073 0.314278 0.181456 0.292683 0.235381 0.226415 0.314523 0.337451 0.2321 0.190692 0.233439 0.30006 0.308556 0.268232 0.212985 0.236908 0.196471 0.362628 0.192465 0.189055 0.176303 0.395448 0.226964 0.243135 0.256469 0.22239 0.253965 0.238164 0.339543 0.276493 0.219117 0.34155 0.214039 0.388673 0.163445 0.284126 0.196305 0.164986 0.31028 0.441618 0.329397 0.256888 0.381206 0.195768 0.318183 0.273093 0.278357 0.169908 0.272998 0.263655 0.254594 0.312809 0.269626 0.067411 0.383223 0.37931 0.191724 0.270489 0.262314 0.239672 0.264436 0.180368 0.216868 0.203412 0.150622 0.290519 0.290775 0.287106 0.349357 0.134866 0.169974 0.278269 0.35454 0.282132 0.286016 0.424946 0.263404 0.222656 0.150417 0.225073 0.193692 0.258708 0.114071 0.130925 0.336195 0.29502 0.280082 0.310482 0.253407 0.157425 0.204624 0.241798 0.274457 0.318812 0.3187 0.372615 0.18969 0.346028 0.159726 0.393855 0.238975 0.268078 0.223484 0.273522 0.225435 0.196356 0.297758 0.196146 0.242638 0.380944 0.271878 0.246516 0.271082 0.272747 0.144801 0.305881 0.244687 0.279726 0.254203 0.349963 0.291943 0.285846 0.156667 0.204597 0.159679 0.322125 0.296537 0.273806 0.268932 0.28281 0.33802 0.275357 0.31787 0.150323 0.199902 0.30576 0.327474 0.27681 0.247274 0.240235 0.178687 0.253333 0.318613 0.323638 0.384694 0.205269 0.19649 0.306002 0.304917 0.23643 0.268539 0.331886 0.236779 0.290295 0.291267 0.300881 0.148994 0.232755 0.182284 0.343466 0.256721 0.137642 0.351024 0.321009 0.265219 0.270182 0.11137 0.266437 0.315657 0.254413 0.306538 0.181177 0.388168 0.291605 0.274684 0.256895 0.154702 0.330561 0.237959 0.268267 0.107554 0.10027 0.348792 0.140657 0.183658 0.409042 0.209212 0.315881 0.260673 0.263849 0.288578 0.28507 0.269495 0.329189 0.232892 0.169297 0.304777 0.237929 0.266842 0.269429 0.137884 0.345301 0.239612 0.304889 0.361275 0.135751 0.32775 0.294365 0.182602 0.23965 0.213716 0.342062 0.130929 0.283312 0.137014 0.291256 0.259146 0.254369 0.281751 0.247639 0.282122 0.288058 0.149106 0.338114 0.220631 0.325427 0.212165 0.216096 0.227356 0.15136 0.219794 0.326398 0.209198 0.347308 0.337329 0.245649 0.158818 0.191251 0.261039 0.333009 0.236431 0.345305 0.125392 0.269505 0.103554 0.212813 0.191758 0.259104 0.35782 0.379457 0.210598 0.293038 0.203985 0.25403 0.138143 0.195161 0.253087 0.303529 0.265479 0.16528 0.407893 0.353072 0.352161 0.076868 0.353085 0.258201 0.217496 0.325697 0.229395 0.250909 0.120578 0.305055 0.247786 0.290669 0.151038 0.113993 0.342005 0.202972 0.373711 0.331318 0.162419 0.137204 0.324989 0.215783 0.248614 0.266795 0.214418 0.169454 0.366715 0.276179 0.196195 0.303732 0.24685 0.101506 0.207445 0.226313 0.336129 0.166082 0.234731 0.231829 0.213456 0.335892 0.276165 0.245573 0.22989 0.328155 0.254986 0.393319 0.381086 0.334885 0.342415 0.222414 0.111957 0.26481 0.146282 0.328967 0.27095 0.315592 0.224611 0.23463 0.225345 0.301244 0.296792 0.254067 0.198736 0.223822 0.333125 0.238111 0.216064 0.276703 0.260344 0.338373 0.19904 0.342957 0.319747 0.231791 0.180249 0.308094 0.285943 0.288524 0.331046 0.207087 0.314249 0.289866 0.399655 0.371902 0.256803 0.3369 0.257767 0.364807 0.3338 0.278142 0.251125 0.265652 0.329067 0.208754 0.147801 0.286244 0.281612 0.168288 0.299315 0.296236 0.308907 0.211772 0.259712 0.092783 0.175842 0.327067 0.266256 0.271321 0.179063 0.344783 0.417265 0.276424 0.188825 0.304554 0.364537 0.196436 0.227253 0.103029 0.090697 0.184244 0.174507 0.285983 0.228055 0.294619 0.152128 0.377603 0.269673 0.24764 0.356992 0.425095 0.24127 0.37376 0.322645 0.271281 0.278966 0.378964 0.19028 0.279396 0.3444 0.230092 0.169752 0.360471 0.323707 0.224373 0.242993 0.272254 0.275359 0.238121 0.234367 0.243334 0.184224 0.20168 0.295325 0.207615 0.150702 0.233537 0.245604 0.159318 0.267932 0.268081 0.343509 0.268002 0.226893 0.178115 0.14906 0.195991 0.324358 0.362353 0.315911 0.28809 0.207899 0.293733 0.240741 0.199769 0.319583 0.240738 0.135804 0.316249 0.186065 0.292504 0.293955 0.241858 0.242963 0.210469 0.242798 0.29241 0.226624 0.282783 0.223408 0.272082 0.298961 0.171414 0.239574 0.322712 0.199328 0.2361 0.285515 0.171895 0.29078 0.320008 0.281675 0.212163 0.225268 0.259012 0.337939 0.311033 0.256802 0.266745 0.338419 0.256853 0.331427 0.268175 0.287825 0.226077 0.31801 0.280507 0.313263 0.20732 0.268962 0.253875 0.323619 0.33861 0.367099 0.133348 0.123619 0.133911 0.233413 0.244742 0.229432 0.361705 0.285955 0.219986 0.365915 0.259044 0.160004 0.167197 0.342843 0.221724 0.283116 0.131451 0.206005 0.244895 0.184127 0.324982 0.296697 0.261492 0.22559 0.213175 0.357466 0.287452 0.396946 0.275889 0.342034 0.329079 0.428746 0.197564 0.22399 0.369657 0.263114 0.17363 0.192614 0.352253 0.212086 0.286075 0.290056 0.32724 0.354505 0.285459 0.226395 0.285548 0.341137 0.320809 0.258471 0.164106 0.268955 0.299149 0.138956 0.248856 0.075799 0.180866 0.23551 0.256479 0.142756 0.292785 0.148926 0.21745 0.192862 0.198022 0.236053 0.348571 0.178439 0.295864 0.202733 0.30248 0.199345 0.27876 0.364695 0.107291 0.395608 0.299366 0.281683 0.333633 0.28298 0.308702 0.274296 0.110046 0.173121 0.134612 0.334517 0.287982 0.281838 0.249734 0.22287 0.261104 0.258176 0.339577 0.239355 0.265431 0.372851 0.374023 0.352915 0.315836 0.297419 0.363123 0.109356 0.349902 0.22437 0.208202 0.204636 0.236493 0.389424 0.222613 0.264403 0.0881 0.281957 0.286607 0.305883 0.279113 0.313179 0.25698 0.345531 0.166 0.264582 0.26718 0.182931 0.327042 0.270755 0.209931 0.363678 0.274955 0.194332 0.323965 0.247372 0.29311 0.344191 0.277835 0.280079 0.308959 0.284846 0.28817 0.315901 0.147761 0.090578 0.282608 0.192041 0.12891 0.208947 0.331414 0.166731 0.274351 0.108481 0.231095 0.255007 0.277631 0.187769 0.161465 0.277225 0.168017 0.333649 0.284211 0.271722 0.198798 0.327833 0.190344 0.098195 0.202239 0.25281 0.379954 0.27357 0.249334 0.250663 0.379407 0.273216 0.161754 0.276386 0.177945 0.09951 0.341649 0.192176 0.232372 0.262084 0.270602 0.170872 0.368869 0.269335 0.314207 0.273101 0.251242 0.294322 0.273105 0.234884 0.229455 0.198928 0.313209 0.241231 0.13036 0.276873 0.258813 0.415687 0.21962 0.232272 0.185559 0.298377 0.223822 0.304961 0.199509 0.280045 0.325789 0.27974 0.345892 0.264825 0.444602 0.215488 0.340069 0.318284 0.234897 0.238441 0.146909 0.29025 0.293415 0.230036 0.252605 0.367212 0.222679 0.380365 0.277449 0.226916 0.154679 0.252383 0.248057 0.285903 0.288777 0.250141 0.193409 0.324723 0.291904 0.119043 0.254168 0.265833 0.349241 0.268284 0.319066 0.323595 0.36481 0.330851 0.262795 0.222544 0.18042 0.247478 0.251899 0.218787 0.171321 0.272558 0.273993 0.322452 0.17707 0.180668 0.144681 0.318053 0.270449 0.185478 0.156946 0.256338 0.276633 0.180233 0.120216 0.239418 0.32611 0.292823 0.231579 0.14875 0.136077 0.141088 0.290758 0.246316 0.302522 0.301468 0.0875 0.344734 0.212703 0.392562 0.34527 0.237906 0.219943 0.226301 0.292899 0.199179 0.221227 0.168732 0.243998 0.258903 0.245216 0.268387 0.206358 0.289803 0.217888 0.226646 0.221378 0.245847 0.304779 0.254241 0.275208 0.27349 0.311238 0.315522 0.27604 0.283399 0.207989 0.270861 0.18186 0.255364 0.357083 0.352821 0.274141 0.301241 0.293817 0.186987 0.267208 0.346323 0.316743 0.317464 0.242138 0.220729 0.286429 0.337721 0.275874 0.32673 0.175294 0.244156 0.218295 0.091239 0.263392 0.182946 0.302183 0.059901 0.317298 0.307821 0.322679 0.189734 0.379938 0.235732 0.291988 0.28035 0.250718 0.22729 0.235256 0.321999 0.204259 0.208679 0.172898 0.202192 0.295275 0.242915 0.144253 0.127714 0.268431 0.265701 0.249934 0.226769 0.180945 0.354696 0.252747 0.287134 0.232 0.405383 0.264349 0.281483 0.174265 0.118297 0.265384 0.323721 0.201041 0.287173 0.15217 0.382841 0.235433 0.316242 0.294598 0.303552 0.398657 0.072892 0.367511 0.26187 0.240952 0.288932 0.381113 0.197241 0.242554 0.246245 0.248554 0.338441 0.289209 0.255694 0.296363 0.220022 0.216313 0.296495 0.383887 0.164658 0.232336 0.201178 0.271929 0.270765 0.268074 0.283902 0.368294 0.224271 0.122846 0.253764 0.192385 0.250519 0.340289 0.379062 0.214486 0.365821 0.185037 0.360796 0.270737 0.307553 0.153913 0.275466 0.214896 0.314369 0.161683 0.219871 0.317078 0.1723 0.229941 0.31448 0.199516 0.304355 0.346169 0.183766 0.275811 0.310974 0.203305 0.263222 0.269103 0.210976 0.353259 0.312411 0.201114 0.147864 0.234655 0.28128 0.306997 0.20008 0.311518 0.189502 0.254927 0.234619 0.209287 0.189025 0.151199 0.30943 0.341948 0.243091 0.340171 0.133583 0.266757 0.247487 0.375916 0.301121 0.24015 0.216313 0.320844 0.272164 0.282852 0.219645 0.09044 0.183534 0.31403 0.283008 0.193589 0.376047 0.12207 0.188805 0.33238 0.250406 0.243237 0.307348 0.179378 0.258408 0.233129 0.342161 0.112658 0.308234 0.294785 0.299547 0.193615 0.174455 0.303865 0.159758 0.133883 0.135245 0.182423 0.285637 0.208162 0.348971 0.316472 0.290929 0.29719 0.394253 0.375194 0.183316 0.313934 0.144721 0.264113 0.309673 0.258356 0.306959 0.203649 0.350442 0.26138 0.21743 0.223952 0.205957 0.16895 0.278309 0.285292 0.285188 0.256936 0.163509 0.296987 0.284596 0.269161 0.127686 0.255464 0.197506 0.234057 0.202981 0.154868 0.245357 0.074878 0.291886 0.251108 0.224155 0.393822 0.160748 0.395014 0.117976 0.333602 0.279401 0.273339 0.311792 0.221833 0.142726 0.366766 0.256705 0.419228 0.234654 0.276844 0.368132 0.303709 0.208796 0.148562 0.179675 0.369331 0.322489 0.275388 0.098124 0.324507 0.265305 0.229702 0.237659 0.38311 0.402087 0.161453 0.337198 0.259134 0.320536 0.34479 0.431261 0.259326 0.204105 0.37503 0.329186 0.241716 0.205905 0.28817 0.19025 0.338149 0.361309 0.315106 0.314905 0.21493 0.176018 0.215224 0.301316 0.278517 0.336381 0.244947 0.308569 0.235248 0.220236 0.342217 0.177676 0.163545 0.099105 0.21781 0.304933 0.205517 0.230638 0.154737 0.380352 0.338214 0.306745 0.210974 0.30772 0.289642 0.321872 0.051023 0.258546 0.25742 0.312499 0.253448 0.195081 0.175637 0.305718 0.261246 0.16628 0.153487 0.216343 0.226766 0.199127 0.246515 0.256892 0.365522 0.332838 0.231304 0.184054 0.201039 0.312976 0.241636 0.305522 0.282438 0.257685 0.381777 0.277787 0.237259 0.197782 0.194587 0.314373 0.157846 0.180315 0.192633 0.32988 0.219774 0.212767 0.337839 0.282007 0.278456 0.303827 0.323001 0.119036 0.16664 0.292284 0.295089 0.173328 0.224361 0.284641 0.259831 0.159575 0.26554 0.24258 0.26502 0.242386 0.36224 0.167873 0.274576 0.232257 0.294756 0.365737 0.229304 0.241702 0.242583 0.315514 0.261597 0.207061 0.291639 0.342618 0.412422 0.212352 0.268583 0.222879 0.347508 0.245471 0.289329 0.264436 0.245213 0.346767 0.392726 0.261484 0.225391 0.236129 0.164651 0.241369 0.210927 0.311418 0.256314 0.343202 0.273142 0.298662 0.275236 0.156157 0.271659 0.276661 0.358709 0.251126 0.367387 0.136207 0.264149 0.2336 0.229585 0.283028 0.182543 0.228378 0.278625 0.303548 0.146436 0.244537 0.294102 0.221011 0.240009 0.086297 0.335798 0.303508 0.24102 0.165346 0.314797 0.238683 0.255459 0.212375 0.251245 0.134348 0.330945 0.253163 0.184662 0.321519 0.237309 0.286613 0.250591 0.300156 0.198211 0.253029 0.168499 0.316019 0.25372 0.337715 0.154117 0.395922 0.201241 0.352337 0.195002 0.268028 0.298354 0.351043 0.365179 0.291613 0.08123 0.218004 0.190302 0.304732 0.1993 0.264798 0.331599 0.247095 0.250834 0.295508 0.227956 0.411904 0.199073 0.175045 0.242089 0.14289 0.301815 0.252956 0.260651 0.284987 0.19406 0.197725 0.139164 0.152633 0.1431 0.247675 0.292625 0.22208 0.169523 0.250312 0.246657 0.260649 0.113637 0.11522 0.349168 0.256983 0.268471 0.198356 0.202347 0.319412 0.221949 0.243146 0.356733 0.107357 0.171368 0.307137 0.096764 0.360548 0.251771 0.21946 0.149717 0.203061 0.238053 0.291521 0.267715 0.273393 0.175583 0.294107 0.209529 0.400832 0.218598 0.292617 0.241775 0.270046 0.250317 0.347085 0.311522 0.410851 0.223259 0.19415 0.241883 0.347414 0.223548 0.25649 0.384433 0.295666 0.230497 0.337077 0.32045 0.314341 0.217293 0.17719 0.288831 0.208642 0.328219 0.343755 0.233994 0.256826 0.361744 0.305028 0.161413 0.177043 0.37553 0.340868 0.28242 0.208614 0.091438 0.259338 0.240517 0.302892 0.208337 0.252851 0.220114 0.274506 0.235173 0.103599 0.400237 0.351032 0.271062 0.252293 0.282518 0.272807 0.288058 0.302551 0.391686 0.233033 0.241363 0.170947 0.175349 0.200568 0.262635 0.200086 0.285427 0.331445 0.250161 0.313744 0.32838 0.274666 0.389737 0.319956 0.233605 0.247686 0.303716 0.281883 0.296829 0.293515 0.195616 0.314866 0.262013 0.321917 0.353672 0.305674 0.261004 0.27738 0.304957 0.30365 0.307659 0.206868 0.280281 0.267933 0.198901 0.230333 0.247269 0.387865 0.298598 0.308174 0.239708 0.176086 0.274485 0.334294 0.255537 0.138845 0.241144 0.157229 0.150327 0.307024 0.259279 0.12711 0.283088 0.23767 0.270749 0.264337 0.246642 0.054888 0.342582 0.242386 0.235189 0.20992 0.309828 0.278824 0.145633 0.236073 0.209734 0.254842 0.280822 0.212673 0.372608 0.387311 0.405267 0.105243 0.133839 0.300548 0.212393 0.236637 0.215289 0.337611 0.245099 0.326612 0.176566 0.256054 0.337925 0.248593 0.315184 0.363329 0.235767 0.18825 0.258562 0.342775 0.265211 0.389259 0.324327 0.192678 0.333222 0.193678 0.189646 0.262241 0.304808 0.15983 0.270437 0.380645 0.277807 0.140092 0.260304 0.230717 0.339922 0.238747 0.24133 0.16255 0.186493 0.330253 0.332985 0.249989 0.241698 0.090905 0.145283 0.241863 0.213543 0.202179 0.29516 0.263584 0.31912 0.265408 0.321702 0.368025 0.196215 0.271665 0.356461 0.302052 0.284014 0.360041 0.209274 0.228221 0.302886 0.278407 0.332992 0.303609 0.115946 0.127638 0.195684 0.325873 0.259683 0.306716 0.2907 0.272883 0.181913 0.199891 0.14 0.290317 0.169167 0.265677 0.326331 0.2543 0.307667 0.342226 0.322612 0.231437 0.402716 0.09149 0.241343 0.139146 0.302752 0.20426 0.334349 0.32319 0.160323 0.306629 0.150082 0.276451 0.256576 0.25176 0.283762 0.22211 0.100308 0.166243 0.184799 0.209278 0.346017 0.322568 0.174785 0.243917 0.364226 0.24326 0.266506 0.287889 0.322826 0.1608 0.350679 0.197102 0.264183 0.17637 0.255667 0.290545 0.265001 0.234506 0.191762 0.361901 0.271669 0.273116 0.303158 0.212263 0.205372 0.294223 0.343432 0.151156 0.314913 0.262313 0.283567 0.095167 0.230677 0.255868 0.329957 0.21485 0.280454 0.347853 0.16378 0.20512 0.297455 0.255915 0.284156 0.286836 0.305902 0.255358 0.339446 0.269981 0.239 0.124803 0.196421 0.224486 0.272238 0.304591 0.10054 0.286416 0.262926 0.280288 0.213935 0.217834 0.317389 0.248542 0.278919 0.314536 0.257652 0.290108 0.179103 0.224682 0.328695 0.319453 0.033888 0.259085 0.284751 0.344655 0.247428 0.105951 0.130812 0.194482 0.285162 0.278633 0.154511 0.243588 0.317473 0.293117 0.183936 0.321006 0.239465 0.239071 0.218927 0.24017 0.27516 0.266365 0.174211 0.321761 0.274727 0.368403 0.23871 0.27362 0.253506 0.221852 0.184912 0.201486 0.292111 0.252431 0.3624 0.362691 0.200528 0.207532 0.335792 0.289141 0.200878 0.24949 0.240226 0.356872 0.293623 0.241195 0.135341 0.238058 0.222526 0.38175 0.196755 0.185845 0.288848 0.178696 0.40602 0.379499 0.163518 0.231013 0.285655 0.285288 0.06933 0.149806 0.29384 0.251007 0.303707 0.215217 0.298774 0.278206 0.253172 0.24193 0.302144 0.379603 0.281801 0.317107 0.323478 0.404051 0.232711 0.267833 0.219642 0.341124 0.299647 0.222882 0.218108 0.23046 0.336851 0.261101 0.379332 0.178262 0.30571 0.193983 0.2966 0.371392 0.185734 0.317554 0.249555 0.341469 0.308546 0.296063 0.230268 0.139944 0.266695 0.321061 0.245233 0.194123 0.26262 0.24945 0.326817 0.36097 0.198479 0.173667 0.291655 0.218113 0.19985 0.287557 0.207053 0.273187 0.263271 0.245527 0.239511 0.298027 0.140379 0.26423 0.102178 0.219486 0.241736 0.353705 0.155 0.349666 0.301896 0.213701 0.170251 0.236404 0.346875 0.142864 0.418117 0.229381 0.236555 0.264732 0.170239 0.320139 0.291449 0.43265 0.333747 0.145538 0.280104 0.273231 0.24845 0.354753 0.280721 0.311082 0.401226 0.251035 0.19408 0.327128 0.308028 0.33973 0.231701 0.091246 0.192009 0.166901 0.332143 0.268012 0.376588 0.379748 0.288462 0.256918 0.118839 0.181191 0.163743 0.334245 0.25053 0.268581 0.356032 0.337952 0.254882 0.269474 0.263291 0.359724 0.214278 0.266507 0.190816 0.290018 0.242863 0.278693 0.162417 0.21734 0.260003 0.298341 0.225568 0.290423 0.212139 0.299654 0.306079 0.317687 0.26045 0.332604 0.24015 0.359621 0.350283 0.295422 0.256032 0.341013 0.25041 0.303307 0.230463 0.228229 0.293623 0.277345 0.196582 0.251226 0.317708 0.344754 0.302374 0.274986 0.199467 0.232957 0.19048 0.189039 0.289108 0.293784 0.318715 0.349082 0.264208 0.215452 0.187465 0.223375 0.309182 0.334306 0.08099 0.369741 0.2858 0.344604 0.393316 0.295897 0.234088 0.255397 0.24502 0.307371 0.289537 0.248621 0.375804 0.271979 0.269855 0.26554 0.317551 0.292574 0.283299 0.16801 0.277242 0.325069 0.278524 0.269174 0.366909 0.307329 0.257109 0.367072 0.247305 0.283166 0.257981 0.246002 0.308111 0.34154 0.208865 0.309396 0.127877 0.355394 0.243447 0.202319 0.267647 0.109216 0.332383 0.318759 0.138421 0.1456 0.38044 0.27045 0.316386 0.2454 0.281826 0.243142 0.246743 0.306522 0.397362 0.24923 0.271763 0.344871 0.318562 0.15374 0.254478 0.27452 0.221004 0.328768 0.165727 0.190939 0.274241 0.153435 0.303504 0.325286 0.327254 0.368225 0.296256 0.221815 0.296162 0.249887 0.284433 0.247804 0.286863 0.24868 0.195943 0.236452 0.166944 0.238105 0.216105 0.228678 0.220119 0.240869 0.252831 0.283414 0.249898 0.351803 0.282664 0.19651 0.129991 0.234863 0.286985 0.355679 0.350887 0.269474 0.342975 0.31639 0.333999 0.26468 0.270898 0.289425 0.376049 0.376233 0.285219 0.176127 0.225401 0.209314 0.188037 0.274954 0.260819 0.229261 0.251537 0.222781 0.297984 0.327384 0.260632 0.217611 0.170765 0.324406 0.206786 0.286147 0.419972 0.314493 0.175023 0.29508 0.333396 0.144056 0.217456 0.469845 0.29842 0.230718 0.180332 0.193768 0.222822 0.362142 0.274108 0.279149 0.203282 0.281102 0.293849 0.290126 0.281559 0.244891 0.345844 0.258129 0.348018 0.355829 0.250224 0.268169 0.296351 0.3056 0.247319 0.219955 0.25604 0.296844 0.304028 0.349995 0.383549 0.293334 0.339749 0.304309 0.324105 0.325343 0.186927 0.271131 0.312428 0.22419 0.179679 0.342639 0.140912 0.233434 0.338969 0.280793 0.314044 0.26706 0.151507 0.387678 0.255027 0.29203 0.224003 0.356749 0.300272 0.212644 0.244642 0.324023 0.149691 0.279488 0.211715 0.203786 0.23213 0.337729 0.217883 0.278563 0.259971 0.281723 0.260064 0.297565 0.291408 0.315763 0.311099 0.17767 0.152164 0.27589 0.182726 0.351927 0.16665 0.400032 0.401201 0.227923 0.242232 0.257438 0.304504 0.257977 0.274556 0.24518 0.273143 0.323112 0.220515 0.283801 0.186811 0.343452 0.292768 0.418744 0.35789 0.197226 0.315788 0.296842 0.195608 0.316673 0.196109 0.324281 0.121133 0.265791 0.256509 0.33277 0.191011 0.149944 0.255311 0.347496 0.19778 0.344703 0.233549 0.18901 0.212371 0.199193 0.256759 0.105796 0.185249 0.278308 0.28235 0.26526 0.204722 0.104485 0.351693 0.115925 0.260836 0.292535 0.352738 0.188577 0.165284 0.219864 0.212597 0.290068 0.262762 0.224452 0.268424 0.28299 0.400604 0.205375 0.128348 0.204622 0.255816 0.294561 0.194067 0.274465 0.224648 0.263058 0.404216 0.231684 0.12605 0.172005 0.29227 0.300684 0.348428 0.35552 0.312214 0.309361 0.140143 0.333925 0.181785 0.220454 0.204207 0.242946 0.273545 0.202497 0.147347 0.302517 0.175256 0.35107 0.167528 0.309888 0.268167 0.21441 0.126328 0.268129 0.226623 0.098671 0.281765 0.247433 0.125041 0.139599 0.305477 0.313639 0.286728 0.311199 0.30224 0.29301 0.238416 0.242349 0.195192 0.107549 0.11232 0.302243 0.270734 0.319873 0.214445 0.291867 0.281876 0.185383 0.228975 0.174474 0.215034 0.200802 0.226338 0.275651 0.274854 0.318228 0.167231 0.295201 0.289411 0.132154 0.269736 0.30073 0.242238 0.207813 0.151522 0.263754 0.323507 0.345658 0.266694 0.150321 0.341103 0.255797 0.209833 0.168533 0.321149 0.335554 0.191251 0.159094 0.125104 0.232847 0.368269 0.379802 0.23277 0.279286 0.21542 0.222321 0.326167 0.260568 0.227752 0.320367 0.313078 0.268257 0.170333 0.327721 0.255482 0.34847 0.264764 0.335961 0.300056 0.327092 0.328593 0.194926 0.385962 0.297753 0.375949 0.380555 0.355899 0.189771 0.21126 0.262277 0.236591 0.201533 0.299954 0.232947 0.329018 0.285766 0.210325 0.185232 0.226463 0.36179 0.33306 0.29936 0.315304 0.338018 0.278036 0.262072 0.155703 0.376961 0.272309 0.253106 0.312292 0.316486 0.341282 0.21312 0.338498 0.322901 0.26892 0.304847 0.309342 0.282444 0.237377 0.303028 0.240979 0.155915 0.120808 0.347 0.401136 0.402679 0.23556 0.324746 0.227112 0.321852 0.376984 0.340947 0.196998 0.205403 0.254134 0.285068 0.21696 0.3139 0.312738 0.23501 0.253062 0.358948 0.233297 0.270095 0.29972 0.283682 0.287656 0.240905 0.215516 0.273289 0.253064 0.144782 0.201743 0.297988 0.332649 0.214471 0.34445 0.268629 0.217099 0.304861 0.097785 0.37173 0.284951 0.105541 0.197369 0.322165 0.285162 0.32428 0.219134 0.299683 0.223979 0.189001 0.258874 0.130247 0.181742 0.123644 0.288737 0.368825 0.294829 0.198983 0.257179 0.205737 0.338997 0.143208 0.320091 0.37076 0.249721 0.209047 0.180504 0.258053 0.386269 0.182066 0.171031 0.27119 0.364741 0.096803 0.328504 0.105523 0.361002 0.183555 0.205098 0.162869 0.143352 0.242475 0.209238 0.238909 0.24488 0.243111 0.222231 0.238165 0.306728 0.303812 0.299125 0.236393 0.237306 0.293108 0.180521 0.116752 0.301126 0.305649 0.313705 0.197909 0.261369 0.231765 0.196678 0.345408 0.264698 0.107518 0.367891 0.331926 0.315325 0.279391 0.198271 0.267594 0.306961 0.259271 0.290005 0.248653 0.287925 0.369392 0.184227 0.125136 0.184272 0.262419 0.248944 0.272544 0.266673 0.169561 0.26839 0.19739 0.227838 0.295442 0.257638 0.374523 0.18419 0.254526 0.224729 0.3614 0.333377 0.147903 0.391539 0.217535 0.344098 0.293959 0.169248 0.229286 0.16898 0.244351 0.313108 0.22198 0.317343 0.358247 0.123024 0.215457 0.162846 0.261714 0.25448 0.30651 0.339513 0.116982 0.275656 0.207976 0.179297 0.343191 0.129958 0.332665 0.255667 0.256064 0.283662 0.319359 0.1732 0.367246 0.310579 0.347084 0.218064 0.315686 0.219428 0.280569 0.150143 0.312953 0.192101 0.221818 0.301581 0.240004 0.246492 0.188191 0.269075 0.355344 0.280336 0.314914 0.309766 0.291004 0.216349 0.207466 0.113665 0.153319 0.348706 0.268561 0.164354 0.266433 0.189125 0.379186 0.206538 0.207418 0.131965 0.335232 0.259156 0.346953 0.262983 0.260644 0.205429 0.366243 0.314088 0.371799 0.317036 0.233571 0.25229 0.221801 0.227257 0.229359 0.271402 0.277665 0.246075 0.244744 0.336343 0.316107 0.299746 0.202515 0.247312 0.279162 0.225854 0.311382 0.330833 0.180685 0.207757 0.370707 0.373258 0.332888 0.321219 0.346532 0.224599 0.327389 0.140054 0.262424 0.098377 0.181016 0.316327 0.224761 0.217186 0.14656 0.174417 0.148514 0.177805 0.285268 0.362398 0.21909 0.272704 0.294552 0.126652 0.232692 0.299226 0.389986 0.335846 0.085405 0.245538 0.223728 0.252878 0.29273 0.130257 0.214848 0.274491 0.253412 0.429266 0.244748 0.260017 0.221666 0.354354 0.330873 0.330763 0.28501 0.259823 0.260463 0.223501 0.322513 0.220808 0.420419 0.26752 0.382789 0.243373 0.268805 0.297081 0.110823 0.154354 0.167202 0.252212 0.197006 0.184236 0.255204 0.289323 0.284083 0.359115 0.285745 0.266233 0.291182 0.318514 0.299165 0.248149 0.295091 0.283382 0.21841 0.362855 0.211378 0.351469 0.190762 0.343211 0.247281 0.22147 0.177843 0.281957 0.373045 0.211688 0.23709 0.30048 0.184875 0.354991 0.353409 0.252608 0.211215 0.31908 0.319692 0.270091 0.154994 0.15979 0.305124 0.283662 0.222604 0.311413 0.352023 0.22601 0.324424 0.176925 0.243407 0.271935 0.353088 0.287249 0.301586 0.3003 0.089891 0.264799 0.254581 0.106587 0.231742 0.194624 0.30515 0.335635 0.192313 0.192068 0.245356 0.207499 0.251873 0.244347 0.185199 0.29971 0.186648 0.140188 0.153983 0.336408 0.284844 0.302465 0.297623 0.293334 0.32748 0.245525 0.314736 0.192732 0.257645 0.25756 0.254779 0.281632 0.129906 0.333073 0.238804 0.329807 0.232916 0.298959 0.283637 0.296343 0.316823 0.234875 0.302889 0.288849 0.290555 0.145605 0.267965 0.23329 0.193781 0.205435 0.278537 0.209011 0.215734 0.297904 0.400461 0.102025 0.410053 0.232733 0.281902 0.275718 0.258549 0.159495 0.301373 0.193821 0.220109 0.285627 0.096634 0.409718 0.160406 0.313771 0.318575 0.336793 0.122787 0.23182 0.232925 0.283046 0.285448 0.099932 0.297065 0.314017 0.309557 0.286464 0.241117 0.223573 0.325126 0.164388 0.191918 0.305319 0.228586 0.197051 0.313327 0.330577 0.324534 0.317898 0.28348 0.190318 0.294405 0.317922 0.32573 0.230377 0.24861 0.170729 0.309124 0.271111 0.286033 0.294271 0.128038 0.219857 0.377212 0.252369 0.0994 0.285456 0.332783 0.332216 0.272299 0.136302 0.201811 0.129009 0.193572 0.257391 0.179745 0.264834 0.369435 0.275661 0.226181 0.255804 0.340098 0.182521 0.3533 0.260662 0.282575 0.21345 0.239423 0.249299 0.155034 0.267733 0.212901 0.236718 0.108874 0.258682 0.117815 0.261753 0.346832 0.251943 0.277696 0.286244 0.297062 0.235485 0.286595 0.336487 0.275009 0.322186 0.188802 0.296897 0.30904 0.189164 0.227726 0.201214 0.17474 0.177556 0.231082 0.3083 0.221449 0.233937 0.300276 0.34783 0.267216 0.201804 0.34376 0.156921 0.183252 0.133241 0.317972 0.200157 0.134736 0.237831 0.313602 0.302371 0.248471 0.353019 0.336938 0.314923 0.301153 0.322639 0.226729 0.333337 0.295949 0.291311 0.254233 0.331688 0.258135 0.401805 0.261532 0.181251 0.346663 0.235692 0.197399 0.225309 0.201397 0.196748 0.224724 0.284864 0.233194 0.21645 0.279745 0.315533 0.199016 0.206695 0.300107 0.227587 0.264288 0.293081 0.295345 0.343057 0.176143 0.201364
)";
    istringstream instring(builtin_distributions);
    ifstream infile;
    
    // read distributions from file if the user has specified a file_path
    if (filepath)
    {
        infile.open(filepath);
        if (!infile.is_open())
        {
            string filepath_str(filepath);
            outError("Error in reading "+filepath_str+": "+ ERR_READ_INPUT);
        }
    }
    
    // select the appropriate source to read distributions
    istream& is = infile.is_open()?static_cast<std::istream&>(infile):instring;
    
    // parse distributions from file's content
    string line;
    
    // parse distributions one by one
    while (getline(is, line))
    {
        // Ignore empty lines
        if (line.length() > 0)
        {
            // validate the input format: <distribution_name> <num_rand_numbers> <random_num_1> <random_num_2> ... <random_num_n>
            int pos = line.find(' ');
            if (pos != std::string::npos)
            {
                // extract distribution_name, and list of random_numbers
                string distribution_name = line.substr(0, pos);
                line.erase(0, pos + 1);
                string random_numbers = line.substr(0, line.length());
                
                // extract num_rand_numbers
                int num_rand_numbers = std::count(random_numbers.begin(), random_numbers.end(), ' ');
                // ignore the space if it's in the last character
                if (random_numbers[random_numbers.length()-1] == ' ')
                    num_rand_numbers--;
                
                // add the new distribution to the list of distributions
                Distribution distribution = {random_numbers, num_rand_numbers + 1};
                Params::getInstance().distributions.insert (std::pair<string,Distribution>(distribution_name, distribution));
            }
            else
                outError("Each distribution should be defined in one line as <distribution_name> <num_rand_numbers> <random_num_1> <random_num_2> ... <random_num_n>");
        }
    }
    
    // close file if it's openning
    if (infile.is_open())
        infile.close();
}

double random_number_from_distribution(string distribution_name, bool non_negative)
{
    // randomly generate a number from a uniform distribution
    if (distribution_name.compare("uniform") == 0)
        return random_double();
        
    Distribution distribution = Params::getInstance().distributions[distribution_name];
    string random_numbers_str = distribution.random_numbers_str;
    
    // check whether distribution_name could be found
    if (random_numbers_str.length() == 0)
    {
        if (distribution_name.length()>0)
            outError("Expecting a double or a distribution name, but found "+distribution_name+". Could not found the distribution named " + distribution_name);
        else
            outError("Expecting a double or a distribution name, but found an empty string");
    }
    
    // attempt up to 1000 times to pick a random (positive) number from the distribution
    double random_number;
    for (int attempt = 0; attempt < 1000; ++attempt)
    {
        // convert random_numbers_str to istringstream
        istringstream iss_random_numbers(random_numbers_str);
        
        // draw a random number
        int rand_index = random_int(distribution.pool_size) + 1;
        
        // extract the selected number from iss_random_numbers
        for (int i = 0; i<rand_index; i++)
            iss_random_numbers >> random_number;
        
        // don't need to retry if the random number meets the output requirement
        if (!non_negative || random_number >= 0)
            break;
    }
    
    // validate the output number
    if (non_negative && random_number < 0)
        outError("Sorry! We failed to generate a random non-negative number from the distribution " + distribution_name + " after 1,000 attempts!");
    
    return random_number;
}

double convert_double_with_distribution(const char *str, bool non_negative)
{
    string input(str);
    // convert the number from the input string if possible
    if (is_number(input))
        return convert_double(str);
    // if a distribution is specified -> randomly generate a number from that distribution
    else
        return random_number_from_distribution(input, non_negative);
}

double convert_double_with_distribution_and_upperbound(string input, double upper_bound, bool non_negative)
{
    double random_double = 0;
    if (is_number(input))
    {
        random_double = convert_double(input.c_str());
        if (random_double >= upper_bound)
            outError("The input number ("+input+") must be less than "+convertDoubleToString(upper_bound)+". Please check and try again!");
    }
    else
        random_double = random_number_from_distribution_with_upperbound(input, upper_bound, non_negative);
    
    return random_double;
}

double random_number_from_distribution_with_upperbound(string distribution_name, double upper_bound, bool non_negative)
{
    double random_double;
    
    // limit 1000 attempts to generate a random number from user-specified distribution
    for (int i = 0; i < 1000; i++)
    {
        random_double = random_number_from_distribution(distribution_name, non_negative);
        if (random_double < upper_bound && random_double >= 0)
            break;
    }
    
    if (random_double >= upper_bound || random_double < 0)
        outError("Unfortunately, could not generate a random number between 0 and "+convertDoubleToString(upper_bound)+" using the list/distribution "+distribution_name+" after 1000 attempts.");
    
    return random_double;
}

void normalize_frequencies(double* freqs, int num_states, double total_freqs, bool show_warning)
{
    ASSERT(num_states > 0);
    // calculate the total_freqs if it's not provided
    if (total_freqs == -1)
    {
        total_freqs = 0;
        for (int i = 0; i < num_states; i++)
            total_freqs += freqs[i];
    }
    
    // normalize the freqs
    if (fabs(total_freqs) < 1e-5)
        outError("Sum of state frequencies must be greater than zero!");
    
    if (fabs(total_freqs-1.0) >= 1e-7)
    {
        total_freqs = 1/total_freqs;
        if (show_warning)
            outWarning("Normalizing state frequencies so that sum of them equals to 1");
        for (int i = 0; i < num_states; i++)
            freqs[i] *= total_freqs;
    }
}

void normalize_frequencies_from_index(double* freqs, int num_states, int starting_index)
{
    ASSERT(num_states > 0);
    // calculate the total_freqs
    double total_freqs = 0;
    for (int i = starting_index; i < starting_index+num_states; i++)
        total_freqs += freqs[i];
    
    // normalize the freqs
    if (fabs(total_freqs) < 1e-5)
        outError("Sum of state frequencies must be greater than zero!");
    total_freqs = 1/total_freqs;
    for (int i = starting_index; i < starting_index+num_states; i++)
        freqs[i] *= total_freqs;
}

bool is_number(const std::string& s)
{
    char* end = nullptr;
    double val = strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0' && val != HUGE_VAL;
}

void random_frequencies_from_distributions(double *freqs, int num_states, string list_distribution_names)
{
    size_t num_comma = std::count(list_distribution_names.begin(), list_distribution_names.end(), ',');
    if (num_comma + 1 != num_states)
        outError("The number of distributions specified in "+list_distribution_names+" ("+convertIntToString(num_comma+1)+") is different from the number of states ("+convertIntToString(num_states)+"). Please check and try again!");
    
    // calculate the total frequency
    double total_freq = 0;
    for (int i = 0; i<num_states; i++)
    {
        // extract the name of distribution used for the current state
        int pos = list_distribution_names.find(',');
        string distribution_name = list_distribution_names.substr(0, pos);
        list_distribution_names.erase(0, pos + 1);
        
        freqs[i] = random_number_from_distribution(distribution_name, true);
        total_freq += freqs[i];
    }
    
    // normalize the frequencies
    for (int i = 0; i<num_states; i++)
        freqs[i] /= total_freq;
}

void readWeightFile(Params &params, int ntaxa, double &scale, StrVector &tax_name, DoubleVector &tax_weight) {
    cout << "Reading scale factor and taxa weights file " << params.param_file << " ..." << endl;
    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(params.param_file);
        string name, tmp;

        in >> tmp;
        scale = convert_double(tmp.c_str());

        for (; !in.eof() && ntaxa > 0; ntaxa--) {
            // remove the failbit
            in.exceptions(ios::badbit);
            if (!(in >> name)) break;
            // set the failbit again
            in.exceptions(ios::failbit | ios::badbit);

            tax_name.push_back(name);
            // read the sequence weight
            in >> tmp;
            tax_weight.push_back(convert_double(tmp.c_str()));
        }
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
    } catch (ios::failure) {
        outError(ERR_READ_INPUT);
    } catch (string str) {
        outError(str);
    }
}

void readStringFile(const char* filename, int max_num, StrVector &strv) {
    try {
        ifstream in;
        // set the failbit and badbit
        in.exceptions(ios::failbit | ios::badbit);
        in.open(filename);
        string name;

        // remove the failbit
        in.exceptions(ios::badbit);
        for (; !in.eof() && max_num > 0; max_num--) {
            if (!(in >> name)) break;
            strv.push_back(name);
        }
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
    } catch (ios::failure) {
        outError(ERR_READ_INPUT);
    }
}

void readInitTaxaFile(Params &params, int ntaxa, StrVector &tax_name) {
    cout << "Reading initial taxa set file " << params.initial_file << " ..." << endl;
    readStringFile(params.initial_file, ntaxa, tax_name);
}

void printString2File(string myString, string filename) {
    ofstream myfile(filename.c_str());
    if (myfile.is_open()) {
        myfile << myString;
        myfile.close();
    } else {
        cout << "Unable to open file " << filename << endl;
    }
}

void readInitAreaFile(Params &params, int nareas, StrVector &area_name) {
    cout << "Reading initial area file " << params.initial_area_file << " ..." << endl;
    readStringFile(params.initial_area_file, nareas, area_name);
}

void readAreasBoundary(char *file_name, MSetsBlock *areas, double *areas_boundary) {

    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(file_name);

        int nset;
        in >> nset;
        if (nset != areas->getNSets())
            throw "File has different number of areas";
        int pos = 0, seq1, seq2;
        for (seq1 = 0; seq1 < nset; seq1++) {
            string seq_name;
            in >> seq_name;
            if (seq_name != areas->getSet(seq1)->name)
                throw "Area name " + seq_name + " is different from " + areas->getSet(seq1)->name;
            for (seq2 = 0; seq2 < nset; seq2++) {
                in >> areas_boundary[pos++];
            }
        }
        // check for symmetric matrix
        for (seq1 = 0; seq1 < nset - 1; seq1++) {
            if (areas_boundary[seq1 * nset + seq1] <= 1e-6)
                throw "Diagonal elements of distance matrix should represent the boundary of single areas";
            for (seq2 = seq1 + 1; seq2 < nset; seq2++)
                if (areas_boundary[seq1 * nset + seq2] != areas_boundary[seq2 * nset + seq1])
                    throw "Shared boundary between " + areas->getSet(seq1)->name + " and " + areas->getSet(seq2)->name + " is not symmetric";
        }


        in.close();
        cout << "Areas relation matrix was read from " << file_name << endl;
    } catch (const char *str) {
        outError(str);
    } catch (string str) {
        outError(str);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, file_name);
    }

}

void readTaxaSets(char *filename, MSetsBlock *sets) {
    TaxaSetNameVector *allsets = sets->getSets();
    try {
        int count = 0;
        ifstream in;
        // set the failbit and badbit
        in.exceptions(ios::failbit | ios::badbit);
        in.open(filename);
        string name;

        // remove the failbit
        in.exceptions(ios::badbit);
        while (!in.eof()) {
            int ntaxa = 0;
            string str;
            if (!(in >> str)) break;
            ntaxa = convert_int(str.c_str());
            if (ntaxa <= 0) throw "Number of taxa must be > 0";
            count++;
            //allsets->resize(allsets->size()+1);
            TaxaSetName *myset = new TaxaSetName;
            allsets->push_back(myset);
            myset->name = "";
            myset->name += count;
            for (; ntaxa > 0; ntaxa--) {
                string str;
                if (!(in >> str)) throw "Cannot read in taxon name";
                if ((ntaxa > 1) && in.eof()) throw "Unexpected end of file while reading taxon names";
                myset->taxlist.push_back(str);
            }
        }
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
        if (count == 0) throw "No set found, you must specify at least 1 set";
    } catch (ios::failure) {
        outError(ERR_READ_INPUT);
    } catch (const char *str) {
        outError(str);
    } catch (string str) {
        outError(str);
    }
}

// Parse the profile mixture model
// MIX{R+Fx} -> MIX{S+FO,S+FO,...,S+FO} with x classes and S is a linked substitution matrix (i.e. linked exchangeabilities)
// OR R+Fx -> MIX{S+FO,S+FO,...,S+FO} with x classes and S is a linked substitution matrix (i.e. linked exchangeabilities)
// and x < 1000
// return true if it is a linked substitution matrix
bool parseProfileMixModelStr(string& model_str) {
    if (model_str.length() == 0)
        return false;

    string modelstr = model_str;
    // change to upper character
    transform(modelstr.begin(), modelstr.end(), modelstr.begin(), ::toupper);

    int pos_mix = modelstr.find("MIX{");
    if (pos_mix != string::npos) {
        int pos_endBrac = modelstr.find_last_of('}');
        if (pos_endBrac == string::npos || pos_endBrac < pos_mix+4) {
            outError(modelstr + " is not in a correct format");
        }
        // get the substring inside MIX{....}
        modelstr = modelstr.substr(pos_mix+4, pos_endBrac-pos_mix-4);
    }

    bool isLinkedSubst = false;
    int pos_F = modelstr.find("+F");
    if (pos_F != string::npos) {
        int endpos = pos_F+2;
        // get the end pos of the integer followed by +F
        while (endpos < modelstr.length() && isdigit(modelstr[endpos]) && modelstr[endpos] != '.')
            endpos++;
        if (endpos >= modelstr.length() || (modelstr[endpos] == '+' || modelstr[endpos] == ',' || modelstr[endpos] == '}')) {
            // Not an integer followed by a character like +F3X4
            if (endpos > pos_F+2 && endpos-pos_F-2 < 4) {
                // +Fx appears, where x is an integer, and x < 1000
                int nclass = atoi(modelstr.substr(pos_F+2,endpos-pos_F-2).c_str());
                if (nclass >= 1) {
                    // this is R+Fx where x is an integer
                    string s_model = modelstr.substr(0, pos_F);
                    string RHAS = "";
                    int pos_plus = modelstr.find("+",pos_F+2);
                    if (pos_plus != string::npos) {
                        RHAS = modelstr.substr(pos_plus);
                    }
                    string mix_model = "";
                    if (nclass > 1) {
                        mix_model.append("MIX{");
                    }
                    for (int i = 0; i < nclass; i++) {
                        if (i > 0)
                            mix_model.append(",");
                        mix_model.append(s_model + "+FO");
                    }
                    if (nclass > 1)
                        mix_model.append("}");
                    mix_model.append(RHAS);
                    model_str = mix_model;
                    isLinkedSubst = true;
                }
            }
        }
    }
    return isLinkedSubst;
}

// Parse the profile mixture model
void parseProfileMix(Params& params) {
    if (parseProfileMixModelStr(params.model_name))
        params.optimize_linked_gtr = true;
    if (parseProfileMixModelStr(params.model_joint))
        params.optimize_linked_gtr = true;
}

void get2RandNumb(const int size, int &first, int &second) {
    // pick a random element
    first = random_int(size);
    // pick a random element from what's left (there is one fewer to choose from)...
    second = random_int(size - 1);
    // ...and adjust second choice to take into account the first choice
    if (second >= first) {
        ++second;
    }
}

void quickStartGuide();

void parseArg(int argc, char *argv[], Params &params) {
    int cnt;
    progress_display::setProgressDisplay(false);
    verbose_mode = VB_MIN;
    params.setDefault();

    // store original params
    for (cnt = 1; cnt < argc; cnt++) {
        params.original_params = params.original_params + argv[cnt] + " ";
    }
    
    for (cnt = 1; cnt < argc; cnt++) {
        try {

            if (strcmp(argv[cnt], "-h") == 0 || strcmp(argv[cnt], "--help") == 0) {
#ifdef IQ_TREE
                usage_iqtree(argv, strcmp(argv[cnt], "--help") == 0);
#else
                usage(argv, false);
#endif
                continue;
            }
            if (strcmp(argv[cnt], "-V") == 0 || strcmp(argv[cnt], "-version") == 0 || strcmp(argv[cnt], "--version") == 0) {
                printCopyright(cout);
                exit(EXIT_SUCCESS);
            }
			if (strcmp(argv[cnt], "-ho") == 0 || strcmp(argv[cnt], "-?") == 0) {
				usage_iqtree(argv, false);
				continue;
			}
			if (strcmp(argv[cnt], "-hh") == 0 || strcmp(argv[cnt], "-hhh") == 0) {
#ifdef IQ_TREE
                usage_iqtree(argv, true);
#else
				usage(argv);
#endif
				continue;
			}
			if (strcmp(argv[cnt], "-v0") == 0) {
				verbose_mode = VB_QUIET;
				continue;
			}
			if (strcmp(argv[cnt], "-v") == 0 || strcmp(argv[cnt], "--verbose") == 0) {
				verbose_mode = VB_MED;
				continue;
			}
			if (strcmp(argv[cnt], "-vv") == 0
					|| strcmp(argv[cnt], "-v2") == 0) {
				verbose_mode = VB_MAX;
				continue;
			}
			if (strcmp(argv[cnt], "-vvv") == 0
					|| strcmp(argv[cnt], "-v3") == 0) {
				verbose_mode = VB_DEBUG;
				continue;
			}
			if (strcmp(argv[cnt], "-k") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -k <num_taxa>";
				convert_range(argv[cnt], params.min_size, params.sub_size,
						params.step_size);
				params.k_representative = params.min_size;
				continue;
			}
			if (strcmp(argv[cnt], "-pre") == 0 || strcmp(argv[cnt], "--prefix") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pre <output_prefix>";
				params.out_prefix = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-pp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pp <pd_proportion>";
				convert_range(argv[cnt], params.min_proportion,
						params.pd_proportion, params.step_proportion);
				if (params.pd_proportion < 0 || params.pd_proportion > 1)
					throw "PD proportion must be between 0 and 1";
				continue;
			}
			if (strcmp(argv[cnt], "-mk") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mk <min_taxa>";
				params.min_size = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bud") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bud <budget>";
				convert_range(argv[cnt], params.min_budget, params.budget,
						params.step_budget);
				continue;
			}
			if (strcmp(argv[cnt], "-mb") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mb <min_budget>";
				params.min_budget = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-o") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -o <taxon>";
				params.root = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-optalg") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -opt_alg <1-BFGS|2-BFGS|EM>";
				params.optimize_alg_freerate = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-optlen") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -optlen <BFGS|EM>";
				params.optimize_alg_mixlen = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "-optalg_gammai") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -optalg_gammai <Brent|BFGS|EM>";
                params.optimize_alg_gammai = argv[cnt];
                continue;
            }
            if (strcmp(argv[cnt], "-optalg_treeweight") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -optalg_treeweight <BFGS|EM>";
                params.optimize_alg_treeweight = argv[cnt];
                continue;
            }

            if (strcmp(argv[cnt], "-optalg_qmix") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -optalg_qmix <BFGS|EM>";
                if(strcmp(argv[cnt], "BFGS") != 0 && strcmp(argv[cnt], "EM") != 0)
                    throw "Invalid option for -optalg_qmix : use 'BFGS' or 'EM'";
                params.optimize_alg_qmix = argv[cnt];
                continue;
            }

            if (strcmp(argv[cnt], "-init_nucl_freq") == 0 || strcmp(argv[cnt], "--init_nucl_freq") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -init_nucl_freq <0|1|2>";
                params.estimate_init_freq = convert_int(argv[cnt]);
                if (params.estimate_init_freq > 2)
                    throw "Use -init_nucl_freq <0|1|2>";
                continue;
            }

            //new options added -JD
            if (strcmp(argv[cnt], "--link-exchange-rates") == 0 || strcmp(argv[cnt], "--link-exchange") == 0) {
                params.optimize_linked_gtr = true;
                params.reset_method = "const";
                params.optimize_alg_qmix = "EM";
                continue;
            }
            if (strcmp(argv[cnt], "--gtr20-model") == 0 || strcmp(argv[cnt], "--init-exchange") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use " + string(argv[cnt-1]) + " <protein model name>";
                params.gtr20_model = argv[cnt];
                transform(params.gtr20_model.begin(), params.gtr20_model.end(), params.gtr20_model.begin(), ::toupper);
                continue;
            }
            if (strcmp(argv[cnt], "--guess-multiplier") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --guess-multiplier <value>";
                params.guess_multiplier = convert_double(argv[cnt]);
                continue;
            }
//            if (strcmp(argv[cnt], "--rates-file") == 0) {
//                params.rates_file = true;
//                continue;
//            }
            if (strcmp(argv[cnt], "--reset-method") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --reset-method <const/random>";
                if(strcmp(argv[cnt], "const") != 0 && strcmp(argv[cnt], "random") != 0)
                    throw "Invalid option for --reset-method : use 'const' or 'random'";
                params.reset_method = argv[cnt];
                continue;
            }

			if (strcmp(argv[cnt], "-root") == 0 || strcmp(argv[cnt], "-rooted") == 0) {
				params.is_rooted = true;
				continue;
			}
            
            if (strcmp(argv[cnt], "--root-dist") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --root-dist <maximum-root-move-distance>";
                params.root_move_dist = convert_int(argv[cnt]);
                continue;
            }

            if (strcmp(argv[cnt], "--root-find") == 0) {
                params.root_find = true;
                continue;
            }

            if (strcmp(argv[cnt], "--root-test") == 0) {
                params.root_test = true;
                continue;
            }
            
			if (strcmp(argv[cnt], "-all") == 0) {
				params.find_all = true;
				continue;
			}
			if (strcmp(argv[cnt], "--greedy") == 0) {
				params.run_mode = RunMode::GREEDY;
				continue;
			}
			if (strcmp(argv[cnt], "-pr") == 0
					|| strcmp(argv[cnt], "--pruning") == 0) {
				params.run_mode = RunMode::PRUNING;
				//continue; } if (strcmp(argv[cnt],"--both") == 0) {
				//params.run_mode = BOTH_ALG;
				continue;
			}
			if (strcmp(argv[cnt], "-e") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -e <file>";
				params.param_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-if") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -if <file>";
				params.initial_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-nni_nr_step") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nni_nr_step <newton_raphson_steps>";
				NNI_MAX_NR_STEP = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-ia") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ia <file>";
				params.initial_area_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-u") == 0) {
				// file containing budget information
				cnt++;
				if (cnt >= argc)
					throw "Use -u <file>";
				params.budget_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-dd") == 0) {
				// compute distribution of PD score on random sets
				cnt++;
				if (cnt >= argc)
					throw "Use -dd <sample_size>";
				params.run_mode = RunMode::PD_DISTRIBUTION;
				params.sample_size = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-ts") == 0) {
				// calculate PD score a taxa set listed in the file
				cnt++;
				//params.run_mode = PD_USER_SET;
				if (cnt >= argc)
					throw "Use -ts <taxa_file>";
				params.pdtaxa_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-bound") == 0) {
				// boundary length of areas
				cnt++;
				if (cnt >= argc)
					throw "Use -bound <file>";
				params.areas_boundary_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-blm") == 0) {
				// boundary length modifier
				cnt++;
				if (cnt >= argc)
					throw "Use -blm <boundary_modifier>";
				params.boundary_modifier = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-dist") == 0
					|| strcmp(argv[cnt], "-d") == 0) {
				// calculate distance matrix from the tree
				params.run_mode = RunMode::CALC_DIST;
				cnt++;
				if (cnt >= argc)
					throw "Use -dist <distance_file>";
				params.dist_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-djc") == 0) {
				params.compute_ml_dist = false;
				continue;
			}
            std::string arg = argv[cnt];
            //Todo; move this up, use == rather than strcmp elsewhere, too.
            if (arg=="-mlnj-only" || arg=="--mlnj-only") {
                params.compute_ml_tree_only = true;
                continue;
            }
			if (strcmp(argv[cnt], "-dobs") == 0) {
				params.compute_obs_dist = true;
				continue;
			}
            if (strcmp(argv[cnt], "-experimental") == 0 || strcmp(argv[cnt], "--experimental") == 0) {
                params.experimental = true;
                continue;
            }
            if (strcmp(argv[cnt], "--no-experimental") == 0) {
                params.experimental = false;
                continue;
            }
			if (strcmp(argv[cnt], "-r") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -r <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = YULE_HARDING;
				continue;
			}
			if (strcmp(argv[cnt], "-rs") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rs <alignment_file>";
				params.tree_gen = YULE_HARDING;
				params.aln_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-rstar") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rstar <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = STAR_TREE;
				continue;
			}
			if (strcmp(argv[cnt], "-ru") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ru <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = UNIFORM;
				continue;
			}
            if (strcmp(argv[cnt], "-rbd") == 0) {
                // get birth_death parameters
                cnt++;
                if (cnt >= argc)
                    throw "Use -rbd {<birth_rate>,<death_rate>} <num_taxa>";
                string bd_params = argv[cnt];
                
                // detect the seperator
                char delimiter = ',';
                if (bd_params.find('/') != std::string::npos)
                    delimiter = '/';
                
                if ((bd_params[0]!='{')
                    || (bd_params[bd_params.length()-1]!='}')
                    || (bd_params.find(delimiter) == std::string::npos))
                    throw "Use -rbd {<birth_rate>,<death_rate>} <num_taxa>";
                params.birth_rate = convert_double(bd_params.substr(1, bd_params.find(delimiter) - 1).c_str());
                if (params.birth_rate <= 0)
                    throw "<birth_rate> must be positive";
                bd_params.erase(0, bd_params.find(delimiter) + 1);
                params.death_rate = convert_double(bd_params.substr(0, bd_params.length()-1).c_str());
                if (params.death_rate < 0 || params.death_rate >= params.birth_rate)
                    throw "<death_rate> must be non-negative and less than <birth_rate>";
                
                // normalize birth_rate & death_rate
                double sum_rate = params.death_rate + params.birth_rate;
                if (fabs(sum_rate-1.0) > 1e-5)
                {
                    outWarning("Normalizing birth rate and death rate so that sum of them is equal to 1.");
                    params.death_rate /= sum_rate;
                    params.birth_rate /= sum_rate;
                }
                
                // get #taxa
                cnt++;
                if (cnt >= argc)
                    throw "Use -rbd {<birth_rate>,<death_rate>} <num_taxa>";
                params.sub_size = convert_int(argv[cnt]);
                params.tree_gen = BIRTH_DEATH;
                continue;
            }
			if (strcmp(argv[cnt], "-rcat") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rcat <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = CATERPILLAR;
				continue;
			}
			if (strcmp(argv[cnt], "-rbal") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rbal <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = BALANCED;
				continue;
			}
            if (strcmp(argv[cnt], "--rand") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -rand UNI | CAT | BAL | NET";
                if (strcmp(argv[cnt], "UNI") == 0)
                    params.tree_gen = UNIFORM;
                else if (strcmp(argv[cnt], "CAT") == 0)
                    params.tree_gen = CATERPILLAR;
                else if (strcmp(argv[cnt], "BAL") == 0)
                    params.tree_gen = BALANCED;
                else if (strcmp(argv[cnt], "NET") == 0)
                    params.tree_gen = CIRCULAR_SPLIT_GRAPH;
                else
                    throw "wrong --rand option";
                continue;
            }
            
            if (strcmp(argv[cnt], "--keep-ident") == 0 || strcmp(argv[cnt], "-keep-ident") == 0) {
                params.ignore_identical_seqs = false;
                continue;
            }
			if (strcmp(argv[cnt], "-rcsg") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rcsg <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = CIRCULAR_SPLIT_GRAPH;
				continue;
			}
			if (strcmp(argv[cnt], "-rpam") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rpam <num_splits>";
				params.num_splits = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-rlen") == 0 || strcmp(argv[cnt], "--rlen") == 0) {
				cnt++;
				if (cnt >= argc - 2)
					throw "Use -rlen <min_len> <mean_len> <max_len>";
				params.min_len = convert_double(argv[cnt]);
				params.mean_len = convert_double(argv[cnt + 1]);
				params.max_len = convert_double(argv[cnt + 2]);
				cnt += 2;
				continue;
			}
			if (strcmp(argv[cnt], "-rzero") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rzero <num_zero_branch>";
				params.num_zero_len = convert_int(argv[cnt]);
				if (params.num_zero_len < 0)
					throw "num_zero_len must not be negative";
				continue;
			}
			if (strcmp(argv[cnt], "-rset") == 0) {
				cnt++;
				if (cnt >= argc - 1)
					throw "Use -rset <overlap> <outfile>";
				params.overlap = convert_int(argv[cnt]);
				cnt++;
				params.pdtaxa_file = argv[cnt];
				params.tree_gen = TAXA_SET;
				continue;
			}
			if (strcmp(argv[cnt], "-rep") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rep <repeated_times>";
				params.repeated_time = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-lim") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lim <pd_limit>";
				params.pd_limit = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-noout") == 0) {
				params.nr_output = 0;
				continue;
			}
			if (strcmp(argv[cnt], "-1out") == 0) {
				params.nr_output = 1;
				continue;
			}
			if (strcmp(argv[cnt], "-oldout") == 0) {
				params.nr_output = 100;
				continue;
			}
			if (strcmp(argv[cnt], "-nexout") == 0) {
				params.nexus_output = true;
				continue;
			}
			if (strcmp(argv[cnt], "-exhaust") == 0) {
				params.run_mode = RunMode::EXHAUSTIVE;
				continue;
			}
			if (strcmp(argv[cnt], "-seed") == 0 || strcmp(argv[cnt], "--seed") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -seed <random_seed>";
				params.ran_seed = abs(convert_int(argv[cnt]));
				continue;
			}
			if (strcmp(argv[cnt], "-pdgain") == 0) {
				params.calc_pdgain = true;
				continue;
			}
			if (strcmp(argv[cnt], "-sup") == 0 || strcmp(argv[cnt], "--support") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sup <target_tree_file>";
				params.second_tree = argv[cnt];
				params.consensus_type = CT_ASSIGN_SUPPORT;
				continue;
			}
			if (strcmp(argv[cnt], "-suptag") == 0 || strcmp(argv[cnt], "--suptag") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -suptag <tagname or ALL>";
				params.support_tag = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "-sup2") == 0) {
                outError("Deprecated -sup2 option, please use --gcf --tree FILE");
            }
            
            if (strcmp(argv[cnt], "--rootstrap") == 0) {
                params.consensus_type = CT_ROOTSTRAP;
                cnt++;
                if (cnt >= argc)
                    throw "Use --rootstrap <user_trees_file>";
                params.treeset_file = argv[cnt];
                continue;
            }
            
            if (strcmp(argv[cnt], "--gcf") == 0) {
                if (params.ancestral_site_concordance != 0)
                    throw "Do not specify both --gcf and --scfl";
				params.consensus_type = CT_ASSIGN_SUPPORT_EXTENDED;
                cnt++;
                if (cnt >= argc)
                    throw "Use --gcf <user_trees_file>";
                params.treeset_file = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "--scf") == 0) {
                if (params.ancestral_site_concordance != 0)
                    throw "Do not specify both --scf and --scfl";
                params.consensus_type = CT_ASSIGN_SUPPORT_EXTENDED;
                cnt++;
                if (cnt >= argc)
                    throw "Use --scf NUM_QUARTETS";
                params.site_concordance = convert_int(argv[cnt]);
                if (params.site_concordance < 1)
                    throw "Positive --scf please";
                continue;
            }
            if (strcmp(argv[cnt], "--ascf") == 0) {
                if (params.consensus_type == CT_ASSIGN_SUPPORT_EXTENDED)
                    throw "Do not specify both --scf and --ascf";
                params.ancestral_site_concordance = 1;
                cnt++;
                if (cnt >= argc)
                    throw "Use --ascf NUM_QUARTETS";
                params.site_concordance = convert_int(argv[cnt]);
                if (params.site_concordance < 1)
                    throw "Positive --ascf please";
                continue;
            }
            if (strcmp(argv[cnt], "--bscf") == 0 || strcmp(argv[cnt], "--scfl") == 0) {
                // UPDATE: sCFL now ignore subtrees with all gaps for a particular site
                if (params.consensus_type == CT_ASSIGN_SUPPORT_EXTENDED)
                    throw "Do not specify --scf or --gcf with --scfl";
                params.ancestral_site_concordance = 2;
                cnt++;
                if (cnt >= argc)
                    throw "Use --scfl NUM_QUARTETS";
                params.site_concordance = convert_int(argv[cnt]);
                if (params.site_concordance < 1)
                    throw "Positive --scfl please";
                continue;
            }
            if (strcmp(argv[cnt], "--scflg") == 0) {
                // OUTDATED: with gaps for historical reason
                if (params.consensus_type == CT_ASSIGN_SUPPORT_EXTENDED)
                    throw "Do not specify --scf or --gcf with --scflg";
                params.ancestral_site_concordance = 3;
                cnt++;
                if (cnt >= argc)
                    throw "Use --scflg NUM_QUARTETS";
                params.site_concordance = convert_int(argv[cnt]);
                if (params.site_concordance < 1)
                    throw "Positive --scflg please";
                continue;
            }
            if (strcmp(argv[cnt], "--scf1") == 0) {
                throw "--scf1 option does not exist. Do you mean --scfl?";
                continue;
            }
            if (strcmp(argv[cnt], "--scf-part") == 0 || strcmp(argv[cnt], "--cf-verbose") == 0) {
                params.site_concordance_partition = true;
                continue;
            }
            if (strcmp(argv[cnt], "--cf-quartet") == 0) {
                params.print_cf_quartets = true;
                continue;
            }
            if (strcmp(argv[cnt], "--df-tree") == 0) {
                params.print_df1_trees = true;
                continue;
            }
            if (strcmp(argv[cnt], "--qic") == 0) {
                params.internode_certainty = 1;
                continue;
            }
			if (strcmp(argv[cnt], "-treew") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -treew <tree_weight_file>";
				params.tree_weight_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-con") == 0 || strcmp(argv[cnt], "--con-tree") == 0) {
				params.consensus_type = CT_CONSENSUS_TREE;
				continue;
			}
			if (strcmp(argv[cnt], "-net") == 0 || strcmp(argv[cnt], "--con-net") == 0) {
				params.consensus_type = CT_CONSENSUS_NETWORK;
                continue;
			}
            
            /**MINH ANH: to serve some statistics on tree*/
			if (strcmp(argv[cnt], "-comp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -comp <treefile>";
				params.consensus_type = COMPARE;
				params.second_tree = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-stats") == 0) {
				params.run_mode = RunMode::STATS;
				continue;
			}
			if (strcmp(argv[cnt], "-gbo") == 0) { //guided bootstrap
				cnt++;
				if (cnt >= argc)
					throw "Use -gbo <site likelihod file>";
				params.siteLL_file = argv[cnt];
				//params.run_mode = GBO;
                continue;
			} // MA
            
			if (strcmp(argv[cnt], "-mprob") == 0) { //compute multinomial distribution probability
				cnt++;
				if (cnt >= argc)
					throw "Use -mprob <ref_alignment>";
				params.second_align = argv[cnt];
				//params.run_mode = MPRO;
                continue;
			} // MA
            
			if (strcmp(argv[cnt], "-min") == 0) {
				params.find_pd_min = true;
				continue;
			}
			if (strcmp(argv[cnt], "-excl") == 0) {
				params.exclusive_pd = true;
				continue;
			}
			if (strcmp(argv[cnt], "-endem") == 0) {
				params.endemic_pd = true;
				continue;
			}
			if (strcmp(argv[cnt], "-compl") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -compl <area_name>";
				params.complement_area = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-cluster") == 0) {
				params.branch_cluster = 4;
				cnt++;
				if (cnt >= argc)
					throw "Use -cluster <taxa_order_file>";
				params.taxa_order_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-taxa") == 0) {
				params.run_mode = RunMode::PRINT_TAXA;
				continue;
			}
			if (strcmp(argv[cnt], "-area") == 0) {
				params.run_mode = RunMode::PRINT_AREA;
				continue;
			}
			if (strcmp(argv[cnt], "-scale") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -scale <scaling_factor>";
				params.scaling_factor = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-scaleg") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -scaleg <gene_scale_factor>";
				params.gene_scale_factor = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-scalebranch") == 0) {
				params.run_mode = RunMode::SCALE_BRANCH_LEN;
				cnt++;
				if (cnt >= argc)
					throw "Use -scalebranch <scaling_factor>";
				params.scaling_factor = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-scalenode") == 0) {
				params.run_mode = RunMode::SCALE_NODE_NAME;
				cnt++;
				if (cnt >= argc)
					throw "Use -scalenode <scaling_factor>";
				params.scaling_factor = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-prec") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -prec <numeric_precision>";
				params.numeric_precision = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-lp") == 0) {
				params.run_mode = RunMode::LINEAR_PROGRAMMING;
				continue;
			}
			if (strcmp(argv[cnt], "-lpbin") == 0) {
				params.run_mode = RunMode::LINEAR_PROGRAMMING;
				params.binary_programming = true;
				continue;
			}
			if (strcmp(argv[cnt], "-qp") == 0) {
				params.gurobi_format = true;
				params.quad_programming = true;
				continue;
			}
			if (strcmp(argv[cnt], "-quiet") == 0 || strcmp(argv[cnt], "--quiet") == 0) {
				verbose_mode = VB_QUIET;
				continue;
			}
			if (strcmp(argv[cnt], "-mult") == 0) {
				params.multi_tree = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bi") == 0 || strcmp(argv[cnt], "--burnin") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bi <burnin_value>";
				params.tree_burnin = convert_int(argv[cnt]);
				if (params.tree_burnin < 0)
					throw "Burnin value must not be negative";
				continue;
			}
			if (strcmp(argv[cnt], "-tm") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -tm <tree_max_count>";
				params.tree_max_count = convert_int(argv[cnt]);
				if (params.tree_max_count < 0)
					throw "tree_max_count must not be negative";
				continue;
			}
			if (strcmp(argv[cnt], "-minsup") == 0 || strcmp(argv[cnt], "--sup-min") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -minsup <split_threshold>";
				params.split_threshold = convert_double(argv[cnt]);
				if (params.split_threshold < 0 || params.split_threshold > 1)
					throw "Split threshold must be between 0 and 1";
				continue;
			}
			if (strcmp(argv[cnt], "-minsupnew") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -minsupnew <split_threshold_1/.../split_threshold_k>";
				params.split_threshold_str = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-tw") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -tw <split_weight_threshold>";
				params.split_weight_threshold = convert_double(argv[cnt]);
				if (params.split_weight_threshold < 0)
					throw "Split weight threshold is negative";
				continue;
			}

            if (strcmp(argv[cnt], "-czb") == 0 || strcmp(argv[cnt], "--polytomy") == 0) {
                params.collapse_zero_branch = true;
                continue;
            }

			if (strcmp(argv[cnt], "-swc") == 0) {
				params.split_weight_summary = SW_COUNT;
				continue;
			}
			if (strcmp(argv[cnt], "-swa") == 0) {
				params.split_weight_summary = SW_AVG_ALL;
				continue;
			}
			if (strcmp(argv[cnt], "-swp") == 0) {
				params.split_weight_summary = SW_AVG_PRESENT;
				continue;
			}
			if (strcmp(argv[cnt], "-iwc") == 0) {
				params.test_input = TEST_WEAKLY_COMPATIBLE;
				continue;
			}
			if (strcmp(argv[cnt], "--aln") == 0 || strcmp(argv[cnt], "--msa") == 0 || strcmp(argv[cnt], "-s") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use --aln, -s <alignment_file>";
				params.aln_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "--sequential") == 0) {
                params.phylip_sequential_format = true;
                continue;
            }
            if (strcmp(argv[cnt], "--symtest") == 0) {
                params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "--bisymtest") == 0) {
                params.symtest = SYMTEST_BINOM;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-only") == 0) {
                params.symtest_only = true;
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-remove-bad") == 0) {
                params.symtest_remove = 1;
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-remove-good") == 0) {
                params.symtest_remove = 2;
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-keep-zero") == 0) {
                params.symtest_keep_zero = true;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-type") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --symtest-type SYM|MAR|INT";
                if (strcmp(argv[cnt], "SYM") == 0)
                    params.symtest_type = 0;
                else if (strcmp(argv[cnt], "MAR") == 0)
                    params.symtest_type = 1;
                else if (strcmp(argv[cnt], "INT") == 0)
                    params.symtest_type = 2;
                else
                    throw "Use --symtest-type SYM|MAR|INT";
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-pval") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --symtest-pval PVALUE_CUTOFF";
                params.symtest_pcutoff = convert_double(argv[cnt]);
                if (params.symtest_pcutoff <= 0 || params.symtest_pcutoff >= 1)
                    throw "--symtest-pval must be between 0 and 1";
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }
            
            if (strcmp(argv[cnt], "--symstat") == 0) {
                params.symtest_stat = true;
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-perm") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --symtest-perm INT";
                params.symtest_shuffle = convert_int(argv[cnt]);
                if (params.symtest_shuffle <= 0)
                    throw "--symtest-perm must be positive";
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "-z") == 0 || strcmp(argv[cnt], "--trees") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -z <user_trees_file>";
				params.treeset_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-zb") == 0 || strcmp(argv[cnt], "--test") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -zb <#replicates>";
				params.topotest_replicates = convert_int(argv[cnt]);
				if (params.topotest_replicates < 1000)
					throw "Please specify at least 1000 replicates";
				continue;
			}
            if (strcmp(argv[cnt], "--estimate-model") == 0) {
                params.topotest_optimize_model = true;
                continue;
            }
			if (strcmp(argv[cnt], "-zw") == 0 || strcmp(argv[cnt], "--test-weight") == 0) {
				params.do_weighted_test = true;
				continue;
			}
			if (strcmp(argv[cnt], "-au") == 0 || strcmp(argv[cnt], "--test-au") == 0) {
				params.do_au_test = true;
				continue;
			}
			if (strcmp(argv[cnt], "-sp") == 0 || strcmp(argv[cnt], "-Q") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sp <partition_file>";
				params.partition_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-spp") == 0 || strcmp(argv[cnt], "-p") == 0 || strcmp(argv[cnt], "--partition") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -p <partition_file>";
				params.partition_file = argv[cnt];
				params.partition_type = BRLEN_SCALE;
                params.opt_gammai = false;
				continue;
			}
			if (strcmp(argv[cnt], "-spj") == 0 || strcmp(argv[cnt], "-q") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -q <partition_file>";
				params.partition_file = argv[cnt];
				params.partition_type = BRLEN_FIX;
                params.optimize_alg_gammai = "Brent";
                params.opt_gammai = false;
				continue;
			}
			if (strcmp(argv[cnt], "-M") == 0) {
                params.partition_type = BRLEN_OPTIMIZE;
                continue;
            }

            if (strcmp(argv[cnt], "-spu") == 0 || strcmp(argv[cnt], "-S") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -spu <partition_file>";
                params.partition_file = argv[cnt];
                params.partition_type = TOPO_UNLINKED;
                params.ignore_identical_seqs = false;
                params.buffer_mem_save = true;
                params.print_splits_nex_file = false;
                continue;
            }
            
            if (strcmp(argv[cnt], "--edge") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --edge equal|scale|unlink";
                if (strcmp(argv[cnt], "equal") == 0)
                    params.partition_type = BRLEN_FIX;
                else if (strcmp(argv[cnt], "scale") == 0)
                    params.partition_type = BRLEN_SCALE;
                else if (strcmp(argv[cnt], "unlink") == 0)
                    params.partition_type = BRLEN_OPTIMIZE;
                else
                    throw "Use --edge equal|scale|unlink";
            }
            
            if (strcmp(argv[cnt], "-rcluster") == 0 || strcmp(argv[cnt], "--rcluster") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rcluster <percent>";
                params.partfinder_rcluster = convert_double(argv[cnt]);
                if (params.partfinder_rcluster < 0 || params.partfinder_rcluster > 100)
                    throw "rcluster percentage must be between 0 and 100";
                params.partition_merge = MERGE_RCLUSTER;
				continue;
            }
            if (strcmp(argv[cnt], "-rclusterf") == 0 || strcmp(argv[cnt], "--rclusterf") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rclusterf <percent>";
                params.partfinder_rcluster = convert_double(argv[cnt]);
                if (params.partfinder_rcluster < 0 || params.partfinder_rcluster > 100)
                    throw "rcluster percentage must be between 0 and 100";
                params.partition_merge = MERGE_RCLUSTERF;
				continue;
            }

            if (strcmp(argv[cnt], "-rcluster-max") == 0 || strcmp(argv[cnt], "--rcluster-max") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rcluster-max <num>";
                params.partfinder_rcluster_max = convert_int(argv[cnt]);
                if (params.partfinder_rcluster_max <= 0)
                    throw "rcluster-max must be between > 0";
                if (params.partfinder_rcluster == 100)
                    params.partfinder_rcluster = 99.9999;
                if (params.partition_merge != MERGE_RCLUSTER && params.partition_merge != MERGE_RCLUSTERF)
                    params.partition_merge = MERGE_RCLUSTERF;
				continue;
            }

            if (strcmp(argv[cnt], "--merge") == 0) {
                if (cnt >= argc-1 || argv[cnt+1][0] == '-') {
                    if (params.partfinder_rcluster == 100)
                        params.partfinder_rcluster = 99.9999;
                    params.partition_merge = MERGE_RCLUSTERF;
                    continue;
                }
                cnt++;
                if (cnt >= argc)
                    throw "Use --merge [none|greedy|rcluster|rclusterf|kmeans]";
                if (strcmp(argv[cnt], "none") == 0)
                    params.partition_merge = MERGE_NONE;
                else if (strcmp(argv[cnt], "greedy") == 0)
                    params.partition_merge = MERGE_GREEDY;
                else if (strcmp(argv[cnt], "rcluster") == 0) {
                    if (params.partfinder_rcluster == 100)
                        params.partfinder_rcluster = 99.9999;
                    params.partition_merge = MERGE_RCLUSTER;
                } else if (strcmp(argv[cnt], "rclusterf") == 0) {
                    if (params.partfinder_rcluster == 100)
                        params.partfinder_rcluster = 99.9999;
                    params.partition_merge = MERGE_RCLUSTERF;
                } else if (strcmp(argv[cnt], "kmeans") == 0)
                    params.partition_merge = MERGE_KMEANS;
                else
                    throw "Use --merge [none|greedy|rcluster|rclusterf|kmeans]";
                continue;
            }

            if (strcmp(argv[cnt], "--merge-model") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --merge-model 1|4|ALL|model1,...,modelK";
                params.merge_models = argv[cnt];
                if (params.partition_merge == MERGE_NONE) {
                    if (params.partfinder_rcluster == 100)
                        params.partfinder_rcluster = 99.9999;
                    params.partition_merge = MERGE_RCLUSTERF;
                    continue;
                }
                continue;
            }

            if (strcmp(argv[cnt], "--merge-rate") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --merge-rate rate1,...,rateK";
                params.merge_rates = argv[cnt];
                if (params.partition_merge == MERGE_NONE) {
                    if (params.partfinder_rcluster == 100)
                        params.partfinder_rcluster = 99.9999;
                    params.partition_merge = MERGE_RCLUSTERF;
                    continue;
                }
                continue;
            }

            if (strcmp(argv[cnt], "--merge-log-rate") == 0) {
                params.partfinder_log_rate = true;
                continue;
            }

            if (strcmp(argv[cnt], "--merge-normal-rate") == 0) {
                params.partfinder_log_rate = false;
                continue;
            }

			if (strcmp(argv[cnt], "-keep_empty_seq") == 0) {
				params.remove_empty_seq = false;
				continue;
			}
			
            if (strcmp(argv[cnt], "--terrace_tphast") == 0) {
#ifdef IQTREE_TERRAPHAST
                params.terrace_analysis_tphast = true;
#else
                    throw "Unsupported command: --terrace_tphast.\n"
                        "Please build IQ-TREE with the USE_TERRAPHAST flag.";
#endif
                continue;
            }

            if (strcmp(argv[cnt], "--terrace") == 0) {
                params.terrace_check = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "--no-terrace") == 0) {
                params.terrace_analysis_tphast = false;
                params.terrace_check = false;
                continue;
            }
            
            if (strcmp(argv[cnt], "-no_terrace") == 0) {
                params.terrace_aware = false;
                params.terrace_analysis_tphast = false;
                params.terrace_check = false;
                continue;
            }
            
            if (strcmp(argv[cnt], "--gentrius") == 0) {
                params.terrace_analysis = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "-g_print") == 0) {
                params.print_terrace_trees = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "-g_print_lim") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -g_print_lim <num_of_trees_to_be_output>";
                params.terrace_print_lim = convert_int(argv[cnt]);
                if(params.terrace_print_lim<=0){
                    throw "Invalid value! Use -g_print_lim <trees_num> with trees_num>0";
                }
                continue;
            }
            
            if (strcmp(argv[cnt], "-g_print_induced") == 0) {
                params.print_induced_trees = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "-g_print_m") == 0) {
                params.print_pr_ab_matrix = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "-g_print_m_o") == 0) {
                params.print_m_overlap = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "-pr_ab_matrix") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -pr_ab_matrix <pr_ab_matrix_file>";
                params.pr_ab_matrix = argv[cnt];
                continue;
            }
            
            if (strcmp(argv[cnt], "-g_query") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -g_query <gentrius_query_set_file>";
                params.terrace_query_set = argv[cnt];
                continue;
            }
            
            if (strcmp(argv[cnt], "-g_rm_leaves") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -g_rm_leaves <number_of_leaves_to_be_inserted>";
                params.terrace_remove_m_leaves = convert_int(argv[cnt]);
                if(params.terrace_remove_m_leaves<=0){
                    throw "Invalid value! Use -g_rm_leaves <leaves_num> with leaves_num>0";
                }
                continue;
            }
            
            
            if (strcmp(argv[cnt], "-g_stop_i") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -g_stop_i <number_of_intermediate_trees_to_stop>";
                params.terrace_stop_intermediate_num = convert_int(argv[cnt]);
                if(params.terrace_stop_intermediate_num<0){
                    throw "Invalid value! Use -g_stop_i <trees_num> with trees_num > 0 to stop after trees_num intermediate trees were generated, or use trees_num == 0 to turn off this stopping rule. Default: 10MLN trees.";
                }
                continue;
            }
            
            if (strcmp(argv[cnt], "-g_stop_t") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -g_stop_t <number_of_stand_trees_to_stop>";
                params.terrace_stop_terrace_trees_num = convert_int(argv[cnt]);
                if(params.terrace_stop_terrace_trees_num<0){
                    throw "Invalid value! Use -g_stop_t <trees_num> with trees_num>0 to stop after trees_num trees from the stand were generated, or use trees_num == 0 to turn off this stopping rule. Default: 1MLN trees.";
                }
                continue;
            }
            
            if (strcmp(argv[cnt], "-g_stop_h") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -g_stop_h <number_of_hours_to_stop>";
                params.terrace_stop_time = convert_double(argv[cnt]);
                if(params.terrace_stop_time<0){
                    throw "Invalid value! Use -g_stop_h <h> with h>0 to stop the run after h hours, or use h == 0 to turn off this stopping rule. Default: 7 days.";
                }
                continue;
            }
            
            if (strcmp(argv[cnt], "-g_non_stop") == 0) {
                params.terrace_non_stop = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "-m_only") == 0) {
                params.matrix_order = true;
                params.print_pr_ab_matrix = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "-gen_all_nni") == 0) {
                params.gen_all_NNI = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "-sf") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sf <ngs_file>";
				params.ngs_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-sm") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sm <ngs_mapped_read_file>";
				params.ngs_mapped_reads = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "--fundi") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --fundi taxa_1,...,taxa_n,proportion";
                string fundi_input = argv[cnt];
                
                // parse input taxon set
                size_t pos = 0;

                // detect the seperator
                char delimiter = ',';
                if (fundi_input.find('/') != std::string::npos)
                    delimiter = '/';
                
                while ((pos = fundi_input.find(delimiter)) != std::string::npos) {
                    params.alisim_fundi_taxon_set.push_back(fundi_input.substr(0, pos));
                    fundi_input.erase(0, pos + 1);
                }
                if (params.alisim_fundi_taxon_set.size() == 0 || fundi_input.length() == 0)
                    throw "Use --fundi taxa_1,...,taxa_n,proportion";
                
                // parse proportion
                if (fundi_input == "estimate") {
                    params.alisim_fundi_proportion = 0.0;
                } else {
                    params.alisim_fundi_proportion = convert_double(fundi_input.c_str());
                    if (params.alisim_fundi_proportion > 1 || params.alisim_fundi_proportion <= 0)
                        throw "Proportion in FunDi model must be positive and not greater than 1";
                }
                
                continue;
            }

            if (strcmp(argv[cnt], "--fundi-init-rho") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --fundi-init-rho <proportion>";
                params.fundi_init_proportion = convert_double(argv[cnt]);
                if (params.fundi_init_proportion >= 1 || params.fundi_init_proportion <= 0)
                    throw "Initial proportion in FunDi model must be positive and smaller than 1";
                continue;
            }

            if (strcmp(argv[cnt], "--fundi-init-branch") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --fundi-init-branch <branch_legth>";
                params.fundi_init_branch_length = convert_double(argv[cnt]);
                if (params.fundi_init_branch_length >= params.max_branch_length || params.fundi_init_branch_length <= 0)
                    throw "Initial branch length in FunDi model must be positive and smaller than 10";
                continue;
            }

            if (strcmp(argv[cnt], "-ngs_gap") == 0) {
				params.ngs_ignore_gaps = false;
				continue;
			}
            if (strcmp(argv[cnt], "--skip-checking-memory") == 0) {
                params.alisim_skip_checking_memory = true;
                continue;
            }
            if (strcmp(argv[cnt], "--no-unaligned") == 0) {
                params.alisim_no_export_sequence_wo_gaps = true;
                continue;
            }
            if (strcmp(argv[cnt], "--sub-level-mixture") == 0) {
                params.alisim_mixture_at_sub_level = true;
                continue;
            }
            if (strcmp(argv[cnt], "--branch-scale") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --branch-scale <SCALE>";
                params.alisim_branch_scale = convert_double(argv[cnt]);
                if (params.alisim_branch_scale <= 0)
                    throw "<SCALE> must be positive!";
                continue;
            }
            if (strcmp(argv[cnt], "--site-rate") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --site-rate MEAN/SAMPLING/MODEL";

                // convert option to uppercase
                string option = argv[cnt];
                transform(option.begin(), option.end(), option.begin(), ::toupper);
                
                if (option.compare("MEAN") == 0)
                    params.alisim_rate_heterogeneity = POSTERIOR_MEAN;
                else if (option.compare("SAMPLING") == 0)
                    params.alisim_rate_heterogeneity = POSTERIOR_DIS;
                else if (option.compare("MODEL") == 0)
                    params.alisim_rate_heterogeneity = UNSPECIFIED;
                else
                    throw "Use --site-rate MEAN/SAMPLING/MODEL";
                continue;
            }
            if (strcmp(argv[cnt], "--site-freq") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --site-freq MEAN/SAMPLING/MODEL";

                // convert option to uppercase
                string option = argv[cnt];
                transform(option.begin(), option.end(), option.begin(), ::toupper);
                
                if (option.compare("MEAN") == 0)
                    params.alisim_stationarity_heterogeneity = POSTERIOR_MEAN;
                else if (option.compare("SAMPLING") == 0)
                    params.alisim_stationarity_heterogeneity = POSTERIOR_DIS;
                else if (option.compare("MODEL") == 0)
                    params.alisim_stationarity_heterogeneity = UNSPECIFIED;
                else
                    throw "Use --site-freq MEAN/SAMPLING/MODEL";
                continue;
            }
            if (strcmp(argv[cnt], "--indel") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --indel <INSERTION_RATE>,<DELETION_RATE>";
                
                string arg = argv[cnt];
                // detect the seperator
                char delimiter = ',';
                if (arg.find('/') != std::string::npos)
                    delimiter = '/';
                
                // validate the input
                size_t pos = arg.find(delimiter);
                if (pos == std::string::npos)
                    throw "Use --indel <INSERTION_RATE>,<DELETION_RATE>";
                
                // get INSERTION_RATIO
                params.alisim_insertion_ratio = convert_double(arg.substr(0, pos).c_str());
                if (params.alisim_insertion_ratio < 0)
                    throw "<INSERTION_RATE> must not be negative.";
                
                // remove "<INSERTION_RATE>,"
                arg.erase(0, pos + 1);
                
                // get DELETION_RATE
                params.alisim_deletion_ratio = convert_double(arg.c_str());
                if (params.alisim_deletion_ratio < 0)
                    throw "<DELETION_RATE> must not be negative.";
                
                continue;
            }
            if (strcmp(argv[cnt], "--indel-size") == 0) {
                cnt++;
                string err_msg = "Use --indel-size <INS_DIS>,<DEL_DIS>. Notes: <INS_DIS>,<DEL_DIS> could be names of user-defined distributions, or GEO{<mean>}, NB{<mean>[/<variance>]}, POW{<double_a>[/<int_max>]}, LAV{<double_a>/<int_max>}, which specifies Geometric, Negative Binomial, Zipfian, and Lavalette distribution, respectively.";
                if (cnt >= argc)
                    throw err_msg;
                
                // seek separator position
                string arg = argv[cnt];
                int pos = -1;
                bool inside_bracket = false;
                for (int i = 0; i < arg.length(); i++)
                {
                    // detect brackets
                    if (arg[i] == '{')
                        inside_bracket = true;
                    else if (arg[i] == '}')
                        inside_bracket = false;
                       
                    // only check separator outside brackets
                    if (!inside_bracket && (arg[i] == ',' || arg[i] == '/'))
                    {
                        pos = i;
                        break;
                    }
                }
                
                // make sure the separator is found
                if (pos == -1)
                    throw err_msg;
                
                // get INSERTION_DISTRIBUTION
                string input = arg.substr(0, pos);
                if (input.length() == 0)
                    throw err_msg;
                params.alisim_insertion_distribution = parseIndelDis(input, "Insertion");
                
                // remove "<INSERTION_DISTRIBUTION>,"
                arg.erase(0, pos + 1);
                
                // get DELETION_DISTRIBUTION
                if (arg.length() == 0)
                    throw err_msg;
                params.alisim_deletion_distribution = parseIndelDis(arg, "Deletion");
                
                continue;
            }
            if (strcmp(argv[cnt], "--write-all") == 0) {
                params.alisim_write_internal_sequences = true;
                continue;
            }
            if (strcmp(argv[cnt], "--only-unroot-tree") == 0) {
                params.alisim_only_unroot_tree = true;
                continue;
            }
            if (strcmp(argv[cnt], "--branch-distribution") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --branch-distribution <distribution_name> to specify a distribution, from which branch lengths will be randomly generated.";
                params.branch_distribution = argv[cnt];
                continue;
            }
            if (strcmp(argv[cnt], "--simulation-thresh") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "--simulation-thresh <threshold>";
                params.alisim_simulation_thresh = convert_double(argv[cnt]);
                if (params.alisim_simulation_thresh < 0 || params.alisim_simulation_thresh > 1)
                    throw "<threshold> must be between 0 and 1. Please check and try again!";
                continue;
            }
			if (strcmp(argv[cnt], "-st") == 0 || strcmp(argv[cnt], "--seqtype") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -st BIN or -st DNA or -st AA or -st CODON or -st MORPH or -st CRXX or -st CFxx.";
                string arg = argv[cnt];
                params.sequence_type = argv[cnt];
                
                // handle MORPH{<#STATE>}
                string ERR_MSG = "Please use MORPH{<#STATE>} to specify the number of states for MORPH. <#STATE> should be between 1 and 32.";
                string t_params = argv[cnt];
                string KEYWORD = "MORPH";
                if ((t_params.length() > KEYWORD.length())
                    && (!t_params.substr(0, KEYWORD.length()).compare(KEYWORD)))
                {
                    // validate the input
                    if ((t_params[KEYWORD.length()]!='{')
                        ||(t_params[t_params.length()-1]!='}'))
                        throw ERR_MSG;
                    
                    // remove "MORPH{"
                    t_params.erase(0, KEYWORD.length() + 1);
                    
                    // remove "}"
                    t_params = t_params.substr(0, t_params.length()-1);
                    
                    // extract num_states
                    params.alisim_num_states_morph = convert_int(t_params.c_str());
                    
                    // validate num_states
                    if (params.alisim_num_states_morph < 1 || params.alisim_num_states_morph > 32)
                        throw ERR_MSG;
                    
                    // set seqtype to MORPH (without {<#STATE>})
                    params.sequence_type = strcpy(new char[KEYWORD.length() + 1], KEYWORD.c_str());
                }
                
                // if (arg.substr(0,2) == "CR") params.pomo_random_sampling = true;
                // if (arg.substr(0,2) == "CF" || arg.substr(0,2) == "CR") {
                //     outWarning("Setting the sampling method and population size with this flag is deprecated.");
                //     outWarning("Please use the model string instead (see `iqtree --help`).");
                //     if (arg.length() > 2) {
                //         int ps = convert_int(arg.substr(2).c_str());
                //         params.pomo_pop_size = ps;
                //         if (((ps != 10) && (ps != 2) && (ps % 2 == 0)) || (ps < 2) || (ps > 19)) {
                //             std::cout << "Please give a correct PoMo sequence type parameter; e.g., `-st CF09`." << std::endl;
                //             outError("Custom virtual population size of PoMo not 2, 10 or any other odd number between 3 and 19.");   
                //         }
                //     }
                // }
				continue;
			}
            
			if (strcmp(argv[cnt], "-starttree") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -starttree BIONJ|PARS|PLLPARS";
                else if (strcmp(argv[cnt], "PARS") == 0)
					params.start_tree = STT_PARSIMONY;
				else if (strcmp(argv[cnt], "PLLPARS") == 0)
					params.start_tree = STT_PLL_PARSIMONY;
                else if (START_TREE_RECOGNIZED(argv[cnt])) {
                    params.start_tree_subtype_name = argv[cnt];
                    params.start_tree = STT_BIONJ;
                }
                else
					throw "Invalid option, please use -starttree with BIONJ or PARS or PLLPARS";
				continue;
			}

			if (strcmp(argv[cnt], "-ao") == 0 || strcmp(argv[cnt], "--out-alignment") == 0 || strcmp(argv[cnt], "--out-aln") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ao <alignment_file>";
				params.aln_output = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-as") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -as <aln_site_list>";
				params.aln_site_list = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-an") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -an <ref_seq_name>";
				params.ref_seq_name = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-af") == 0 || strcmp(argv[cnt], "--out-format") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -af phy|fasta";
                string format = argv[cnt];
                transform(format.begin(), format.end(), format.begin(), ::toupper);
				if (strcmp(format.c_str(), "PHY") == 0)
					params.aln_output_format = IN_PHYLIP;
				else if (strcmp(format.c_str(), "FASTA") == 0)
					params.aln_output_format = IN_FASTA;
                else if (strcmp(format.c_str(), "NEXUS") == 0)
                    params.aln_output_format = IN_NEXUS;
                else if (strcmp(format.c_str(), "MAPLE") == 0)
                    params.aln_output_format = IN_MAPLE;
				else
					throw "Unknown output format";
				continue;
			}
            if (strcmp(argv[cnt], "-pathogen") == 0 || strcmp(argv[cnt], "--pathogen") == 0) {
                params.inference_alg = ALG_AUTO;
                continue;
            }
            if (strcmp(argv[cnt], "-pathogen-force") == 0 || strcmp(argv[cnt], "--pathogen-force") == 0) {
                params.inference_alg = ALG_CMAPLE;
                continue;
            }
            if (strcmp(argv[cnt], "--sprta") == 0 ||
                strcmp(argv[cnt], "-sprta") == 0) {
              params.compute_SPRTA = true;

              continue;
            }
            if (strcmp(argv[cnt], "--sprta-zero-branch") == 0 ||
                strcmp(argv[cnt], "-sprta-zero-branch") == 0) {
              params.SPRTA_zero_branches = true;

              continue;
            }
            if (strcmp(argv[cnt], "--sprta-other-places") == 0 ||
                strcmp(argv[cnt], "-sprta-other-places") == 0) {
              params.out_alter_spr = true;

              continue;
            }
            if (strcmp(argv[cnt], "--out-csv") == 0) {
                params.output_format = FORMAT_CSV;
                continue;
            }
            
            if (strcmp(argv[cnt], "--out-tsv") == 0) {
                params.output_format = FORMAT_TSV;
                continue;
            }            

            if (strcmp(argv[cnt], "--figtree") == 0) {
                params.newick_extended_format = true;
                continue;
            }

            if (strcmp(argv[cnt], "-am") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -am <gap_masked_aln>";
				params.gap_masked_aln = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-ac") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ac <concatenate_aln>";
				params.concatenate_aln = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-nogap") == 0) {
				params.aln_nogaps = true;
				continue;
			}
			if (strcmp(argv[cnt], "-noconst") == 0) {
				params.aln_no_const_sites = true;
				continue;
			}
			if (strcmp(argv[cnt], "-alninfo") == 0) {
				params.print_aln_info = true;
				continue;
			}
//			if (strcmp(argv[cnt], "-parstree") == 0) {
				// maximum parsimony
//				params.parsimony_tree = true;
//            continue; } if (strcmp(argv[cnt], "-pars") == 0) {
//                // maximum parsimony
//                params.parsimony = true;
//				continue;
//			}
			if (strcmp(argv[cnt], "-spr") == 0) {
				// subtree pruning and regrafting
				params.tree_spr = true;
				continue;
			}
			if (strcmp(argv[cnt], "-krep") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -krep <num_k>";
				params.k_representative = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-pdel") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pdel <probability>";
				params.p_delete = convert_double(argv[cnt]);
				if (params.p_delete < 0.0 || params.p_delete > 1.0)
					throw "Probability of deleting a leaf must be between 0 and 1";
				continue;
			}
			if (strcmp(argv[cnt], "-pers") == 0 || strcmp(argv[cnt], "--perturb") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pers <perturbation_strength>";
				params.initPS = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-n") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -n <#iterations>";
                if (params.gbo_replicates != 0) {
                    throw("Ultrafast bootstrap does not work with -n option");
                }
				params.min_iterations = convert_int(argv[cnt]);
				params.stop_condition = SC_FIXED_ITERATION;
//                params.autostop = false;
				continue;
			}
			if (strcmp(argv[cnt], "-nparam") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nparam <#iterations>";
				params.num_param_iterations = convert_int(argv[cnt]);
				if (params.num_param_iterations < 0)
					throw "Number of parameter optimization iterations (-nparam) must be non negative";
				continue;
			}

			if (strcmp(argv[cnt], "-nb") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nb <#bootstrap_replicates>";
				params.min_iterations = convert_int(argv[cnt]);
				params.iqp_assess_quartet = IQP_BOOTSTRAP;
//				params.avoid_duplicated_trees = true;
				continue;
			}
			if (strcmp(argv[cnt], "--model") == 0 || strcmp(argv[cnt], "-m") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use --model <model_name>";
                
                params.model_name = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "--length-ratio") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --length-ratio <LENGTH_RATIO>";
                params.alisim_length_ratio = convert_double(argv[cnt]);
                if (params.alisim_length_ratio <= 1)
                    throw "<LENGTH_RATIO> must be greater than 1.";
                continue;
            }
            if (strcmp(argv[cnt], "--init-model") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --init-model FILE";
                params.model_name_init = argv[cnt];
                continue;
            }
            if (strcmp(argv[cnt], "--loop-model") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --loop-model NUM";
                params.model_opt_steps = convert_int(argv[cnt]);
                continue;
            }
            if (strcmp(argv[cnt], "--nonrev-model") == 0 || strcmp(argv[cnt], "-nonrev-model") == 0) {
                params.contain_nonrev = true;
                continue;
            }
			if (strcmp(argv[cnt], "-mset") == 0 || strcmp(argv[cnt], "--mset") == 0 || strcmp(argv[cnt], "--models") == 0 || strcmp(argv[cnt], "-mexchange") == 0 || strcmp(argv[cnt], "--mexchange") == 0 ) {
				cnt++;
				if (cnt >= argc)
					throw "Use " + string(argv[cnt-1]) + " <model_set>";
				params.model_set = argv[cnt];
                if (params.model_set == "non-reversible")
                    params.contain_nonrev = true;
				continue;
			}
			if (strcmp(argv[cnt], "-madd") == 0 || strcmp(argv[cnt], "--madd") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -madd <extra_model_set>";
				params.model_extra_set = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-msub") == 0 || strcmp(argv[cnt], "--msub") == 0 || strcmp(argv[cnt], "--model-sub") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -msub <model_subset>";
				params.model_subset = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-mfreq") == 0 || strcmp(argv[cnt], "--freqs") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mfreq <state_freq_set>";
				params.state_freq_set = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-mrate") == 0 || strcmp(argv[cnt], "--mrate") == 0 || strcmp(argv[cnt], "--rates") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mrate <rate_set>";
				params.ratehet_set = argv[cnt];
				continue;
			}
            
            if (strcmp(argv[cnt], "--score-diff") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --score-diff <score>";
                if (iEquals(argv[cnt], "all"))
                    params.score_diff_thres = -1.0;
                else
                    params.score_diff_thres = convert_double(argv[cnt]);
                continue;
            }
            
			if (strcmp(argv[cnt], "-mdef") == 0 || strcmp(argv[cnt], "--mdef") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mdef <model_definition_file>";
				params.model_def_file = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "--modelomatic") == 0) {
                params.modelomatic = true;
                continue;
            }
			if (strcmp(argv[cnt], "-mredo") == 0 || strcmp(argv[cnt], "--mredo") == 0 || strcmp(argv[cnt], "--model-redo") == 0) {
				params.model_test_again = true;
				continue;
			}
			if (strcmp(argv[cnt], "-mtree") == 0 || strcmp(argv[cnt], "--mtree") == 0) {
				params.model_test_and_tree = 1;
				continue;
			}
			if (strcmp(argv[cnt], "-mretree") == 0) {
				params.model_test_and_tree = 2;
				continue;
			}
			if (strcmp(argv[cnt], "-msep") == 0) {
				params.model_test_separate_rate = true;
				continue;
			}
			if (strcmp(argv[cnt], "-mwopt") == 0 || strcmp(argv[cnt], "--mix-opt") == 0) {
				params.optimize_mixmodel_weight = true;
				continue;
			}
			if (strcmp(argv[cnt], "-mfopt") == 0 || strcmp(argv[cnt], "--mfopt") == 0) {
				params.optimize_mixmodel_freq = true;
				continue;
			}
			if (strcmp(argv[cnt], "--opt-rate-mat") == 0) {
				params.optimize_rate_matrix = true;
				continue;
			}
            if (strcmp(argv[cnt], "-parallel-over-sites") == 0 || strcmp(argv[cnt], "--parallel-over-sites") == 0) {
                params.parallel_over_sites = true;
                continue;
            }

            // parallelization ordered by threads
            if (strcmp(argv[cnt], "-parallel-order-thread") == 0 || strcmp(argv[cnt], "--parallel-order-thread") == 0) {
                params.order_by_threads = true;
                continue;
            }



//			if (strcmp(argv[cnt], "-mh") == 0) {
//				params.mvh_site_rate = true;
//				params.discard_saturated_site = false;
//				params.SSE = LK_NORMAL;
//				continue;
//			}
//			if (strcmp(argv[cnt], "-mhs") == 0) {
//				params.mvh_site_rate = true;
//				params.discard_saturated_site = true;
//				params.SSE = LK_NORMAL;
//				continue;
//			}
			if (strcmp(argv[cnt], "-rl") == 0) {
				params.rate_mh_type = false;
				continue;
			}
			if (strcmp(argv[cnt], "-nr") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nr <mean_rate>";
				params.mean_rate = convert_double(argv[cnt]);
				if (params.mean_rate < 0)
					throw "Wrong mean rate for MH model";
				continue;
			}
			if (strcmp(argv[cnt], "-mstore") == 0) {
				params.store_trans_matrix = true;
				continue;
			}
			if (strcmp(argv[cnt], "-nni_lh") == 0) {
				params.nni_lh = true;
				continue;
			}
			if (strcmp(argv[cnt], "-lmd") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lmd <lambda>";
                params.lambda = convert_double(argv[cnt]);
				if (params.lambda > 1.0)
					throw "Lambda must be in (0,1]";
				continue;
			}

            if (strcmp(argv[cnt], "-lk") == 0) {
				cnt++;
				if (cnt >= argc)
                    throw "-lk x86|SSE|AVX|FMA|AVX512";
                if (strcmp(argv[cnt], "x86") == 0)
                    params.SSE = LK_386;
                else if (strcmp(argv[cnt], "SSE") == 0)
                    params.SSE = LK_SSE2;
                else if (strcmp(argv[cnt], "AVX") == 0)
                    params.SSE = LK_AVX;
                else if (strcmp(argv[cnt], "FMA") == 0)
                    params.SSE = LK_AVX_FMA;
                else if (strcmp(argv[cnt], "AVX512") == 0)
                    params.SSE = LK_AVX512;
                else
                    throw "Incorrect -lk likelihood kernel option";
				continue;
			}

			if (strcmp(argv[cnt], "-safe") == 0 || strcmp(argv[cnt], "--safe") == 0) {
				params.lk_safe_scaling = true;
				continue;
			}

			if (strcmp(argv[cnt], "-safe-seq") == 0) {
				cnt++;
				if (cnt >= argc)
                    throw "-safe-seq <number of sequences>";
				params.numseq_safe_scaling = convert_int(argv[cnt]);
                if (params.numseq_safe_scaling < 10)
                    throw "Too small -safe-seq";
				continue;
			}

            if (strcmp(argv[cnt], "--kernel-nonrev") == 0) {
                params.kernel_nonrev = true;
                continue;
            }

			if (strcmp(argv[cnt], "-f") == 0) {
				cnt++;
				if (cnt >= argc)
				        throw "Use -f <c | o | u | q | ry | ws | mk | <digits>>";
				if (strcmp(argv[cnt], "q") == 0 || strcmp(argv[cnt], "EQ") == 0)
					params.freq_type = FREQ_EQUAL;
				else if (strcmp(argv[cnt], "c") == 0
						|| strcmp(argv[cnt], "EM") == 0)
					params.freq_type = FREQ_EMPIRICAL;
				else if (strcmp(argv[cnt], "o") == 0
						|| strcmp(argv[cnt], "ES") == 0)
					params.freq_type = FREQ_ESTIMATE;
				else if (strcmp(argv[cnt], "u") == 0
						|| strcmp(argv[cnt], "UD") == 0)
					params.freq_type = FREQ_USER_DEFINED;
				else if (strcmp(argv[cnt], "ry") == 0
						|| strcmp(argv[cnt], "RY") == 0)
					params.freq_type = FREQ_DNA_RY;
				else if (strcmp(argv[cnt], "ws") == 0
						|| strcmp(argv[cnt], "WS") == 0)
					params.freq_type = FREQ_DNA_WS;
				else if (strcmp(argv[cnt], "mk") == 0
						|| strcmp(argv[cnt], "MK") == 0)
					params.freq_type = FREQ_DNA_MK;
				else
				        // throws error message if can't parse
				        params.freq_type = parseStateFreqDigits(argv[cnt]);
				continue;
			}

            if (strcmp(argv[cnt], "--keep-zero-freq") == 0) {
                params.keep_zero_freq = true;
                continue;
            }

            if (strcmp(argv[cnt], "--min-freq") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --min-freq NUM";
                params.min_state_freq = convert_double(argv[cnt]);
                if (params.min_state_freq <= 0)
                    throw "--min-freq must be positive";
                if (params.min_state_freq >= 1.0)
                    throw "--min-freq must be < 1.0";
                continue;
            }

            if (strcmp(argv[cnt], "--inc-zero-freq") == 0) {
                params.keep_zero_freq = false;
                continue;
            }

			if (strcmp(argv[cnt], "-fs") == 0 || strcmp(argv[cnt], "--site-freq") == 0) {
                if (params.tree_freq_file)
                    throw "Specifying both -fs and -ft not allowed";
				cnt++;
				if (cnt >= argc)
					throw "Use -fs <site_freq_file>";
				params.site_freq_file = argv[cnt];
//				params.SSE = LK_EIGEN;
				continue;
			}
			if (strcmp(argv[cnt], "-ft") == 0 || strcmp(argv[cnt], "--tree-freq") == 0) {
                if (params.site_freq_file)
                    throw "Specifying both -fs and -ft not allowed";
                cnt++;
				if (cnt >= argc)
					throw "Use -ft <treefile_to_infer_site_frequency_model>";
                params.tree_freq_file = argv[cnt];
                if (params.print_site_state_freq == WSF_NONE)
                    params.print_site_state_freq = WSF_POSTERIOR_MEAN;
                continue;
            }

			if (strcmp(argv[cnt], "-fconst") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -fconst <const_pattern_frequencies>";
				params.freq_const_patterns = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "--nrate") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -c <#rate_category>";
				params.num_rate_cats = convert_int(argv[cnt]);
				if (params.num_rate_cats < 1)
					throw "Wrong number of rate categories";
				continue;
			}
			if (strcmp(argv[cnt], "-cmin") == 0 || strcmp(argv[cnt], "--cmin") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -cmin <#min_rate_category>";
				params.min_rate_cats = convert_int(argv[cnt]);
				if (params.min_rate_cats < 2)
					throw "Wrong number of rate categories for -cmin";
				continue;
			}
			if (strcmp(argv[cnt], "-cmax") == 0 || strcmp(argv[cnt], "--cmax") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -cmax <#max_rate_category>";
				params.max_rate_cats = convert_int(argv[cnt]);
				if (params.max_rate_cats < 2)
					throw "Wrong number of rate categories for -cmax";
				continue;
			}
            if (strcmp(argv[cnt], "-start_subst") == 0 || strcmp(argv[cnt], "--start_subst") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -start_subst <substitution matrix + freq>";
                params.start_subst = argv[cnt];
                continue;
            }
            if (strcmp(argv[cnt], "-skip-opt-combin-subst") == 0 || strcmp(argv[cnt], "--skip-opt-combin-subst") == 0) {
                params.check_combin_q_mat = false;
                continue;
            }
            if (strcmp(argv[cnt], "-qmax") == 0 || strcmp(argv[cnt], "--qmax") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -qmax <#max_mix_classes>";
                params.max_mix_cats = convert_int(argv[cnt]);
                if (params.max_mix_cats < 2)
                    throw "Wrong number of classes in mixture for -qmax. Must be at least 2";
                continue;
            }
            if (strcmp(argv[cnt], "-mrate-twice") == 0 || strcmp(argv[cnt], "--mrate-twice") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -mrate-twice <0|1>";
                int in_option = convert_int(argv[cnt]);
                if (in_option < 0 || in_option > 1)
                    throw "Wrong option for -mrate-twice. Only 0 or 1 is allowed.";
                if (in_option == 1)
                    params.opt_rhas_again = true;
                continue;
            }
            if (strcmp(argv[cnt], "-lrt") == 0 || strcmp(argv[cnt], "--lrt") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -lrt <p-value threshold>";
                params.opt_qmix_pthres = convert_double(argv[cnt]);
                if (params.opt_qmix_pthres <= 0.0 || params.opt_qmix_pthres > 1.0) {
                    throw "Wrong p-value threshold for -opt_qmix_pthres. Must be between 0.0 and 1.0";
                } else {
                    params.opt_qmix_criteria = 1;
                }
                /*
                if (params.opt_qmix_pthres == 0)
                    params.opt_qmix_criteria = 2; // using information critera instead of likelihood-ratio test for estimation of number of classes for Q-Mixture model
                */
                continue;
            }
			if (strcmp(argv[cnt], "-a") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -a <gamma_shape>";
				params.gamma_shape = convert_double(argv[cnt]);
				if (params.gamma_shape <= 0)
					throw "Wrong gamma shape parameter (alpha)";
				continue;
			}

			if (strcmp(argv[cnt], "-amin") == 0 || strcmp(argv[cnt], "--alpha-min") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -amin <min_gamma_shape>";
				params.min_gamma_shape = convert_double(argv[cnt]);
				if (params.min_gamma_shape <= 0)
					throw "Wrong minimum gamma shape parameter (alpha)";
				continue;
			}

			if (strcmp(argv[cnt], "-gmean") == 0 || strcmp(argv[cnt], "--gamma-mean") == 0) {
				params.gamma_median = false;
				continue;
			}
			if (strcmp(argv[cnt], "-gmedian") == 0 || strcmp(argv[cnt], "--gamma-median") == 0) {
				params.gamma_median = true;
				continue;
			}
			if (strcmp(argv[cnt], "-i") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -i <p_invar_sites>";
				params.p_invar_sites = convert_double(argv[cnt]);
				if (params.p_invar_sites < 0)
					throw "Wrong number of proportion of invariable sites";
				continue;
			}
			if (strcmp(argv[cnt], "-optfromgiven") == 0) {
				params.optimize_from_given_params = true;
				continue;
			}
            if (strcmp(argv[cnt], "-hmmster") == 0) {
                params.optimize_params_use_hmm = true;
                params.optimize_params_use_hmm_sm = true;
                params.optimize_params_use_hmm_gm = false;
                params.optimize_params_use_hmm_tm = false;
                params.treemix_optimize_methods = "hmm";
                continue;
            }
            if (strcmp(argv[cnt], "-hmmster{sm}") == 0) {
                params.optimize_params_use_hmm = true;
                params.optimize_params_use_hmm_sm = true;
                params.optimize_params_use_hmm_gm = false;
                params.optimize_params_use_hmm_tm = false;
                params.treemix_optimize_methods = "hmm";
                continue;
            }
            if (strcmp(argv[cnt], "-hmmster{gm}") == 0) {
                params.optimize_params_use_hmm = true;
                params.optimize_params_use_hmm_sm = false;
                params.optimize_params_use_hmm_gm = true;
                params.optimize_params_use_hmm_tm = false;
                params.treemix_optimize_methods = "hmm";
                continue;
            }
            if (strcmp(argv[cnt], "-hmmster{tm}") == 0) {
                params.optimize_params_use_hmm = true;
                params.optimize_params_use_hmm_sm = false;
                params.optimize_params_use_hmm_gm = false;
                params.optimize_params_use_hmm_tm = true;
                params.treemix_optimize_methods = "hmm";
                continue;
            }
//            if (strcmp(argv[cnt], "-hmmonly") == 0) {
//                params.proceed_MAST_after_HMMSTER = false;
//                continue;
//            }
            if (strcmp(argv[cnt], "-hmm_min_stran") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -hmm_min_stran <minimum HMM same-category transition probability>";
                params.HMM_min_stran = convert_double(argv[cnt]);
                if (params.HMM_min_stran >= 1.0 || params.HMM_min_stran < 0.0)
                    throw "Wrong probability for -hmm_min_stran";
                continue;
            }
            if (strcmp(argv[cnt], "-hmm_no_avg_brlen") == 0) {
                params.HMM_no_avg_brlen = true;
                continue;
            }
            if (strcmp(argv[cnt], "-tmix_opt_method") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -tmix_opt_method <hmm/hmm2mast/mast2hmm/mast>";
                params.treemix_optimize_methods = argv[cnt];
                if (strcmp(argv[cnt], "hmm") != 0 && strcmp(argv[cnt], "hmm2mast") != 0 && strcmp(argv[cnt], "mast") != 0 && strcmp(argv[cnt], "mast2hmm") != 0 && strcmp(argv[cnt], "hmast") != 0)
                    throw "Wrong value for -tmix_opt_method";
                continue;
            }
			if (strcmp(argv[cnt], "-brent") == 0) {
				params.optimize_by_newton = false;
				continue;
			}
			if (strcmp(argv[cnt], "-jointopt") == 0) {
				params.optimize_model_rate_joint = true;
				continue;
			}
			if (strcmp(argv[cnt], "-brent_ginvar") == 0) {
				params.optimize_model_rate_joint = false;
				continue;
			}
			if (strcmp(argv[cnt], "-fixbr") == 0 || strcmp(argv[cnt], "-blfix") == 0) {
				params.fixed_branch_length = BRLEN_FIX;
                params.optimize_alg_gammai = "Brent";
                params.opt_gammai = false;
                params.min_iterations = 0;
                params.stop_condition = SC_FIXED_ITERATION;
				continue;
			}
			if (strcmp(argv[cnt], "-blscale") == 0) {
				params.fixed_branch_length = BRLEN_SCALE;
                params.optimize_alg_gammai = "Brent";
                params.opt_gammai = false;
                params.min_iterations = 0;
                params.stop_condition = SC_FIXED_ITERATION;
				continue;
			}
			if (strcmp(argv[cnt], "-blmin") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -blmin <min_branch_length>";
				params.min_branch_length = convert_double(argv[cnt]);
				if (params.min_branch_length < 0.0)
					throw("Negative -blmin not allowed!");
				if (params.min_branch_length == 0.0)
					throw("Zero -blmin is not allowed due to numerical problems");
				if (params.min_branch_length > 0.1)
					throw("-blmin must be < 0.1");

				continue;
			}
			if (strcmp(argv[cnt], "-blmax") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -blmax <max_branch_length>";
				params.max_branch_length = convert_double(argv[cnt]);
				if (params.max_branch_length < 0.5)
					throw("-blmax smaller than 0.5 is not allowed");
				continue;
			}
            if (strcmp(argv[cnt], "--show-lh") == 0) {
                params.ignore_identical_seqs = false;
                params.fixed_branch_length = BRLEN_FIX;
                params.optimize_alg_gammai = "Brent";
                params.opt_gammai = false;
                params.min_iterations = 0;
                params.stop_condition = SC_FIXED_ITERATION;
                verbose_mode = VB_DEBUG;
                params.ignore_checkpoint = true;
                continue;
            }
			if (strcmp(argv[cnt], "-sr") == 0) {
				params.stop_condition = SC_WEIBULL;
				cnt++;
				if (cnt >= argc)
					throw "Use -sr <#max_iteration>";
				params.max_iterations = convert_int(argv[cnt]);
				if (params.max_iterations <= params.min_iterations)
					throw "Specified max iteration must be greater than min iteration";
				continue;
			}
			if (strcmp(argv[cnt], "-nm") == 0 || strcmp(argv[cnt], "--nmax") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nm <#max_iteration>";
				params.max_iterations = convert_int(argv[cnt]);
				if (params.max_iterations <= params.min_iterations)
					throw "Specified max iteration must be greater than min iteration";
				continue;
			}
			if (strcmp(argv[cnt], "-sc") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sc <stop_confidence_value>";
				params.stop_confidence = convert_double(argv[cnt]);
				if (params.stop_confidence <= 0.5
						|| params.stop_confidence >= 1)
					throw "Stop confidence value must be in range (0.5,1)";
				continue;
			}
            if (strcmp(argv[cnt], "--runs") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --runs <number_of_runs>";
                params.num_runs = convert_int(argv[cnt]);
                if (params.num_runs < 1)
                    throw "Positive --runs please";
                continue;
            }
			if (strcmp(argv[cnt], "-gurobi") == 0) {
				params.gurobi_format = true;
				continue;
			}
			if (strcmp(argv[cnt], "-gthreads") == 0) {
				params.gurobi_format = true;
				cnt++;
				if (cnt >= argc)
					throw "Use -gthreads <gurobi_threads>";
				params.gurobi_threads = convert_int(argv[cnt]);
				if (params.gurobi_threads < 1)
					throw "Wrong number of threads";
				continue;
			}
			if (strcmp(argv[cnt], "-b") == 0 || strcmp(argv[cnt], "--boot") == 0 ||
                strcmp(argv[cnt], "-j") == 0 || strcmp(argv[cnt], "--jack") == 0 ||
                strcmp(argv[cnt], "-bo") == 0 || strcmp(argv[cnt], "--bonly") == 0) {
				params.multi_tree = true;
				if (strcmp(argv[cnt], "-bo") == 0 || strcmp(argv[cnt], "--bonly") == 0)
					params.compute_ml_tree = false;
				else
					params.consensus_type = CT_CONSENSUS_TREE;
                if ((strcmp(argv[cnt], "-j") == 0 || strcmp(argv[cnt], "--jack") == 0) && params.jackknife_prop == 0.0)
                    params.jackknife_prop = 0.5;
				cnt++;
				if (cnt >= argc)
					throw "Use -b <num_bootstrap_samples>";
				params.num_bootstrap_samples = convert_int(argv[cnt]);
				if (params.num_bootstrap_samples < 1)
					throw "Wrong number of bootstrap samples";
				if (params.num_bootstrap_samples == 1)
					params.compute_ml_tree = false;
				if (params.num_bootstrap_samples == 1)
					params.consensus_type = CT_NONE;
				continue;
			}
			if (strcmp(argv[cnt], "--bsam") == 0 || strcmp(argv[cnt], "-bsam") == 0 || strcmp(argv[cnt], "--sampling") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bsam <bootstrap_specification>";
				params.bootstrap_spec = argv[cnt];
                params.remove_empty_seq = false;
				continue;
			}
            
            if (strcmp(argv[cnt], "--subsample") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --subsample NUM";
                params.subsampling = convert_int(argv[cnt]);
                continue;
            }
            
            if (strcmp(argv[cnt], "--subsample-seed") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --subsample-seed <random_seed>";
                params.subsampling_seed = convert_int(argv[cnt]);
                continue;
            }

#ifdef USE_BOOSTER
            if (strcmp(argv[cnt], "--tbe") == 0) {
                params.transfer_bootstrap = 1;
                continue;
            }

            if (strcmp(argv[cnt], "--tbe-raw") == 0) {
                params.transfer_bootstrap = 2;
                continue;
            }
#endif

            if (strcmp(argv[cnt], "-bc") == 0 || strcmp(argv[cnt], "--bcon") == 0) {
				params.multi_tree = true;
				params.compute_ml_tree = false;
				cnt++;
				if (cnt >= argc)
					throw "Use -bc <num_bootstrap_samples>";
				params.num_bootstrap_samples = convert_int(argv[cnt]);
				if (params.num_bootstrap_samples < 1)
					throw "Wrong number of bootstrap samples";
				if (params.num_bootstrap_samples > 1)
					params.consensus_type = CT_CONSENSUS_TREE;
				continue;
			}
			if (strcmp(argv[cnt], "-iqppars") == 0) {
				params.iqp_assess_quartet = IQP_PARSIMONY;
				continue;
			}
			if (strcmp(argv[cnt], "-iqp") == 0) {
				params.iqp = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wct") == 0) {
				params.write_candidate_trees = true;
				continue;
			}

			if (strcmp(argv[cnt], "-wt") == 0 || strcmp(argv[cnt], "--treels") == 0) {
				params.write_intermediate_trees = 1;
				continue;
			}

            if (strcmp(argv[cnt], "-wdt") == 0) {
                params.writeDistImdTrees = true;
                continue;
            }

            if (strcmp(argv[cnt], "-wtc") == 0) {
                params.write_intermediate_trees = 1;
                params.print_tree_lh = true;
                continue;
            }

			if (strcmp(argv[cnt], "-wt2") == 0) {
				params.write_intermediate_trees = 2;
//				params.avoid_duplicated_trees = true;
				params.print_tree_lh = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wt3") == 0) {
				params.write_intermediate_trees = 3;
//				params.avoid_duplicated_trees = true;
				params.print_tree_lh = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wbl") == 0) {
				params.print_branch_lengths = true;
				continue;
			}
            if (strcmp(argv[cnt], "-wit") == 0) {
                params.write_init_tree = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "--write-branches") == 0) {
                params.write_branches = true;
                continue;
            }
            
//			if (strcmp(argv[cnt], "-nodup") == 0) {
//				params.avoid_duplicated_trees = true;
//				continue;
//			}
			if (strcmp(argv[cnt], "-rf_all") == 0 || strcmp(argv[cnt], "--tree-dist-all") == 0) {
				params.rf_dist_mode = RF_ALL_PAIR;
				continue;
			}
			if (strcmp(argv[cnt], "-rf_adj") == 0) {
				params.rf_dist_mode = RF_ADJACENT_PAIR;
				continue;
			}
			if (strcmp(argv[cnt], "-rf") == 0 || strcmp(argv[cnt], "--tree-dist") == 0) {
				params.rf_dist_mode = RF_TWO_TREE_SETS;
				cnt++;
				if (cnt >= argc)
					throw "Use -rf <second_tree>";
				params.second_tree = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "-rf1") == 0 || strcmp(argv[cnt], "--tree-dist1") == 0) {
                params.rf_dist_mode = RF_TWO_TREE_SETS;
                params.rf_same_pair = true;
                cnt++;
                if (cnt >= argc)
                    throw "Use --tree-dist1 <second_tree>";
                params.second_tree = argv[cnt];
                continue;
            }
			if (strcmp(argv[cnt], "-rf2") == 0 || strcmp(argv[cnt], "--tree-dist2") == 0) {
				params.rf_dist_mode = RF_TWO_TREE_SETS_EXTENDED;
				cnt++;
				if (cnt >= argc)
					throw "Use -rf2 <second_tree>";
				params.second_tree = argv[cnt];
				continue;
			}
            
            if (strcmp(argv[cnt], "--normalize-dist") == 0) {
                params.normalize_tree_dist = true;
                continue;
            }
            
			if (strcmp(argv[cnt], "-aLRT") == 0) {
				cnt++;
				if (cnt + 1 >= argc)
					throw "Use -aLRT <threshold%> <#replicates>";
				params.aLRT_threshold = convert_int(argv[cnt]);
				if (params.aLRT_threshold < 85 || params.aLRT_threshold > 101)
					throw "aLRT threshold must be between 85 and 100";
				cnt++;
				params.aLRT_replicates = convert_int(argv[cnt]);
				if (params.aLRT_replicates < 1000
						&& params.aLRT_replicates != 0)
					throw "aLRT replicates must be at least 1000";
				continue;
			}
			if (strcmp(argv[cnt], "-alrt") == 0 || strcmp(argv[cnt], "--alrt") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -alrt <#replicates | 0>";
                int reps = convert_int(argv[cnt]);
                if (reps == 0)
                    params.aLRT_test = true;
                else {
                    params.aLRT_replicates = reps;
                    if (params.aLRT_replicates < 1000)
                        throw "aLRT replicates must be at least 1000";
                }
				continue;
			}
			if (strcmp(argv[cnt], "-abayes") == 0 || strcmp(argv[cnt], "--abayes") == 0) {
				params.aBayes_test = true;
				continue;
			}
			if (strcmp(argv[cnt], "-lbp") == 0 || strcmp(argv[cnt], "--lbp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lbp <#replicates>";
				params.localbp_replicates = convert_int(argv[cnt]);
				if (params.localbp_replicates < 1000
						&& params.localbp_replicates != 0)
					throw "Local bootstrap (LBP) replicates must be at least 1000";
				continue;
			}
			if (strcmp(argv[cnt], "-wsl") == 0 || strcmp(argv[cnt], "--sitelh") == 0) {
				params.print_site_lh = WSL_SITE;
				continue;
			}

			if (strcmp(argv[cnt], "-wpl") == 0 || strcmp(argv[cnt], "--partlh") == 0) {
				params.print_partition_lh = true;
				continue;
			}

            if (strcmp(argv[cnt], "-wmp") == 0) {
                params.print_marginal_prob = true;
                continue;
            }

			if (strcmp(argv[cnt], "-wslg") == 0 || strcmp(argv[cnt], "-wslr") == 0) {
				params.print_site_lh = WSL_RATECAT;
				continue;
			}

			if (strcmp(argv[cnt], "-wslm") == 0) {
				params.print_site_lh = WSL_MIXTURE;
				continue;
			}
			if (strcmp(argv[cnt], "-wslmr") == 0 || strcmp(argv[cnt], "-wslrm") == 0) {
				params.print_site_lh = WSL_MIXTURE_RATECAT;
				continue;
			}

			if (strcmp(argv[cnt], "-wspr") == 0) {
				params.print_site_prob = WSL_RATECAT;
				continue;
			}

			if (strcmp(argv[cnt], "-wspm") == 0) {
				params.print_site_prob = WSL_MIXTURE;
				continue;
			}
			if (strcmp(argv[cnt], "-wspmr") == 0 || strcmp(argv[cnt], "-wsprm") == 0) {
				params.print_site_prob = WSL_MIXTURE_RATECAT;
				continue;
			}

			if (strcmp(argv[cnt], "-asr") == 0 || strcmp(argv[cnt], "--ancestral") == 0) {
				params.print_ancestral_sequence = AST_MARGINAL;
                params.ignore_identical_seqs = false;
				continue;
			}

			if (strcmp(argv[cnt], "-asr-min") == 0 || strcmp(argv[cnt], "--asr-min") == 0) {
                cnt++;
				if (cnt >= argc)
					throw "Use -asr-min <probability>";
                
                params.min_ancestral_prob = convert_double(argv[cnt]);
                if (params.min_ancestral_prob < 0 || params.min_ancestral_prob > 1)
                    throw "Minimum ancestral probability [-asr-min] must be between 0 and 1.0";
                continue;
            }

			if (strcmp(argv[cnt], "-asr-joint") == 0) {
				params.print_ancestral_sequence = AST_JOINT;
                params.ignore_identical_seqs = false;
				continue;
			}

			if (strcmp(argv[cnt], "-wsr") == 0 || strcmp(argv[cnt], "--rate") == 0) {
				params.print_site_rate |= 1;
				continue;
			}

            if (strcmp(argv[cnt], "--mlrate") == 0) {
                params.print_site_rate |= 2;
                continue;
            }

            if (strcmp(argv[cnt], "-wsptrees") == 0) {
				params.print_trees_site_posterior = 1;
				continue;
			}
			if (strcmp(argv[cnt], "-wsf") == 0) {
				params.print_site_state_freq = WSF_POSTERIOR_MEAN;
				continue;
			}
			if (strcmp(argv[cnt], "--freq-max") == 0 || strcmp(argv[cnt], "-fmax") == 0) {
				params.print_site_state_freq = WSF_POSTERIOR_MAX;
				continue;
			}
			if (strcmp(argv[cnt], "-wba") == 0) {
				params.print_bootaln = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wbsf") == 0) {
				params.print_boot_site_freq = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wsa") == 0) {
				params.print_subaln = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wtl") == 0) {
				params.print_tree_lh = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wpi") == 0) {
				params.print_partition_info = true;
				params.print_conaln = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wca") == 0) {
				params.print_conaln = true;
				continue;
			}

			if (strcmp(argv[cnt], "-wsplits") == 0) {
				params.print_splits_file = true;
				continue;
			}
            if (strcmp(argv[cnt], "--no-splits.nex") == 0) {
                params.print_splits_nex_file = false;
                continue;
            }
			if (strcmp(argv[cnt], "-ns") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ns <num_simulations>";
				params.whtest_simulations = convert_int(argv[cnt]);
				if (params.whtest_simulations < 1)
					throw "Wrong number of simulations for WH-test";
				continue;
			}
			if (strcmp(argv[cnt], "-mr") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mr <rate_file>";
				params.rate_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-cat_mean") == 0) {
				params.mcat_type |= MCAT_MEAN;
				continue;
			}
			if (strcmp(argv[cnt], "-cat_nolog") == 0) {
				params.mcat_type &= (127 - MCAT_LOG);
				continue;
			}
			if (strcmp(argv[cnt], "-cat_site") == 0) {
				params.mcat_type &= (127 - MCAT_PATTERN);
				continue;
			}
			if (strcmp(argv[cnt], "-tina") == 0) {
				params.do_pars_multistate = true;
                params.ignore_checkpoint = true;
				continue;
			}
			if (strcmp(argv[cnt], "-pval") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pval <gene_pvalue_file>";
				params.gene_pvalue_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-nnitest") == 0) {
				params.testNNI = true;
				continue;
			}
			if (strcmp(argv[cnt], "-anni") == 0) {
				params.approximate_nni = true;
				continue;
			}
			if (strcmp(argv[cnt], "-nnicut") == 0) {
				params.estimate_nni_cutoff = true;
				//nni_cutoff = -5.41/2;
				continue;
			}
			if (strcmp(argv[cnt], "-nnichi2") == 0) {
				params.nni_cutoff = -5.41 / 2;
				continue;
			}
			if (strcmp(argv[cnt], "-nnicutval") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nnicutval <log_diff_value>";
				params.nni_cutoff = convert_double(argv[cnt]);
				if (params.nni_cutoff >= 0)
					throw "cutoff value for -nnicutval must be negative";
				continue;
			}
			if (strcmp(argv[cnt], "-nnisort") == 0) {
				params.nni_sort = true;
				continue;
			}
			if (strcmp(argv[cnt], "-plog") == 0) {
				params.gene_pvalue_loga = true;
				continue;
			}
			if (strcmp(argv[cnt], "-dmp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -dmp <ncbi_taxid>";
				params.ncbi_taxid = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-dmplevel") == 0
					|| strcmp(argv[cnt], "-dmprank") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -dmprank <ncbi_taxon_rank>";
				params.ncbi_taxon_level = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-dmpignore") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -dmpignore <ncbi_ignore_level>";
				params.ncbi_ignore_level = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-dmpname") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -dmpname <ncbi_names_file>";
				params.ncbi_names_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-eco") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -eco <eco_dag_file>";
				params.eco_dag_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-k%") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -k% <k in %>";
				//convert_range(argv[cnt], params.k_percent, params.sub_size, params.step_size);
				params.k_percent = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-diet") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -diet <d in %>";
				convert_range(argv[cnt], params.diet_min, params.diet_max,
						params.diet_step);
				//params.diet = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-up") == 0) {
				params.upper_bound = true;
				continue;
			}
			if (strcmp(argv[cnt], "-upNNI") == 0) {
 				params.upper_bound_NNI = true;
			}
			if (strcmp(argv[cnt], "-upFrac") == 0) {
				cnt++;
				if (cnt >= argc)
				  throw "Use -upFrac <fraction>";
				params.upper_bound_frac = convert_double(argv[cnt]);
			}
			if (strcmp(argv[cnt], "-ecoR") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ecoR <run number>";
				params.eco_run = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bb") == 0 || strcmp(argv[cnt], "-B") == 0 || strcmp(argv[cnt], "--ufboot") == 0 ||
                strcmp(argv[cnt], "-J") == 0 || strcmp(argv[cnt], "--ufjack") == 0) {
                if ((strcmp(argv[cnt], "-J") == 0 || strcmp(argv[cnt], "--ufjack") == 0) && params.jackknife_prop == 0.0)
                    params.jackknife_prop = 0.5;
				cnt++;
				if (cnt >= argc)
					throw "Use -B <#replicates>";
                if (params.stop_condition == SC_FIXED_ITERATION) {
                    throw("Ultrafast bootstrap does not work with -fast, -te or -n option");
                }
				params.gbo_replicates = convert_int(argv[cnt]);
//				params.avoid_duplicated_trees = true;
				if (params.gbo_replicates < 1000)
					throw "#replicates must be >= 1000";
				params.consensus_type = CT_CONSENSUS_TREE;
				params.stop_condition = SC_BOOTSTRAP_CORRELATION;
				//params.nni5Branches = true;
				continue;
			}
			if (strcmp(argv[cnt], "-beps") == 0 || strcmp(argv[cnt], "--beps") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -beps <epsilon>";
				params.ufboot_epsilon = convert_double(argv[cnt]);
				if (params.ufboot_epsilon <= 0.0)
					throw "Epsilon must be positive";
				continue;
			}
			if (strcmp(argv[cnt], "-wbt") == 0 || strcmp(argv[cnt], "--wbt") == 0 || strcmp(argv[cnt], "--boot-trees") == 0) {
				params.print_ufboot_trees = 1;
				continue;
			}
			if (strcmp(argv[cnt], "-wbtl") == 0 || strcmp(argv[cnt], "--wbtl") == 0) {
                // print ufboot trees with branch lengths
				params.print_ufboot_trees = 2;
				continue;
			}
			if (strcmp(argv[cnt], "-bs") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bs <begin_sampling_size>";
				params.check_gbo_sample_size = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bmax") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bmax <max_candidate_trees>";
				params.max_candidate_trees = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bcor") == 0 | strcmp(argv[cnt], "--bcor") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bcor <min_correlation>";
				params.min_correlation = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "--bnni") == 0 || strcmp(argv[cnt], "-bnni") == 0) {
				params.ufboot2corr = true;
                // print ufboot trees with branch lengths
//				params.print_ufboot_trees = 2; // Diep: relocate to be below this for loop
				continue;
			}
			if (strcmp(argv[cnt], "-u2c_nni5") == 0) {
				params.u2c_nni5 = true;
				continue;
			}

			if (strcmp(argv[cnt], "-nstep") == 0 || strcmp(argv[cnt], "--nstep") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nstep <step_iterations>";
				params.step_iterations = convert_int(argv[cnt]);
				if (params.step_iterations < 10
						|| params.step_iterations % 2 == 1)
					throw "At least step size of 10 and even number please";
				params.min_iterations = params.step_iterations;
				continue;
			}
			if (strcmp(argv[cnt], "-boff") == 0) {
				params.online_bootstrap = false;
				continue;
			}
//			if (strcmp(argv[cnt], "-nostore") == 0
//					|| strcmp(argv[cnt], "-memsave") == 0) {
//				params.store_candidate_trees = false;
//				continue;
//			}
            
            if (strcmp(argv[cnt], "--jack-prop") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --jack-prop jackknife_proportion";
                params.jackknife_prop = convert_double(argv[cnt]);
                if (params.jackknife_prop <= 0.0 || params.jackknife_prop >= 1.0)
                    throw "Jackknife proportion must be between 0.0 and 1.0";
                continue;
            }
            
            if (strcmp(argv[cnt], "--robust-phy") == 0) {
                if (params.robust_median)
                    throw "Can't couple --robust-phy with --robust-median";
                cnt++;
                if (cnt >= argc)
                    throw "Use --robust-phy proportion_of_best_sites_to_keep";
                params.robust_phy_keep = convert_double(argv[cnt]);
                if (params.robust_phy_keep <= 0.0 || params.robust_phy_keep > 1.0)
                    throw "--robust-phy parameter must be between 0 and 1";
                // TODO: use Brent (instead of Newton) optimisation of branch lengths
                params.optimize_by_newton = false;
                params.optimize_alg_gammai = "Brent";
                params.optimize_alg_freerate = "2-BFGS";
                continue;
            }

            if (strcmp(argv[cnt], "--robust-median") == 0) {
                if (params.robust_phy_keep < 1.0)
                    throw "Can't couple --robust-phy with --robust-median";
                params.robust_median = true;
                // TODO: use Brent (instead of Newton) optimisation of branch lengths
                params.optimize_by_newton = false;
                params.optimize_alg_gammai = "Brent";
                params.optimize_alg_freerate = "2-BFGS";
                continue;
            }

			if (strcmp(argv[cnt], "-mem") == 0 || strcmp(argv[cnt], "--mem") == 0) {
				cnt++;
				if (cnt >= argc)
                    throw "Use -mem max_mem_size";
				params.lh_mem_save = LM_MEM_SAVE;
                int end_pos;
                double mem = convert_double(argv[cnt], end_pos);
                if (mem < 0)
                    throw "-mem must be non-negative";
                if (argv[cnt][end_pos] == 'G') {
                    params.max_mem_size = mem * 1073741824.0;
                } else if (argv[cnt][end_pos] == 'M') {
                    params.max_mem_size = mem * 1048576.0;
                } else if (argv[cnt][end_pos] == '%'){
                    params.max_mem_size = mem * 0.01;
                    if (params.max_mem_size > 1)
                        throw "-mem percentage must be between 0 and 100";
                } else {
                    if (mem > 1)
                        throw "Invalid -mem option. Example: -mem 200M, -mem 10G";
                    params.max_mem_size = mem;
                }
				continue;
			}
            if (strcmp(argv[cnt], "--save-mem-buffer") == 0) {
                params.buffer_mem_save = true;
                continue;
            }
            if (strcmp(argv[cnt], "--no-save-mem-buffer") == 0) {
                params.buffer_mem_save = false;
                continue;
            }
//			if (strcmp(argv[cnt], "-storetrees") == 0) {
//				params.store_candidate_trees = true;
//				continue;
//			}
			if (strcmp(argv[cnt], "-nodiff") == 0) {
				params.distinct_trees = false;
				continue;
			}
			if (strcmp(argv[cnt], "-treediff") == 0) {
				params.distinct_trees = true;
				continue;
			}
			if (strcmp(argv[cnt], "-norell") == 0) {
				params.use_rell_method = false;
				continue;
			}
			if (strcmp(argv[cnt], "-elw") == 0) {
				params.use_elw_method = true;
				continue;
			}
			if (strcmp(argv[cnt], "-noweight") == 0) {
				params.use_weighted_bootstrap = false;
				continue;
			}
			if (strcmp(argv[cnt], "-nomore") == 0) {
				params.use_max_tree_per_bootstrap = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bweight") == 0) {
				params.use_weighted_bootstrap = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bmore") == 0) {
				params.use_max_tree_per_bootstrap = false;
				continue;
			}
			if (strcmp(argv[cnt], "-gz") == 0) {
				params.do_compression = true;
				continue;
			}
			if (strcmp(argv[cnt], "-newheu") == 0) {
				params.new_heuristic = true;
				// Enable RAxML kernel
				continue;
			}
			if (strcmp(argv[cnt], "-maxtime") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -maxtime <time_in_minutes>";
				params.maxtime = convert_double(argv[cnt]);
				params.min_iterations = 1000000;
				params.stop_condition = SC_REAL_TIME;
				continue;
			}
			if (strcmp(argv[cnt], "--ninit") == 0 || strcmp(argv[cnt], "-ninit") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ninit <number_of_parsimony_trees>";
				params.numInitTrees = convert_int(argv[cnt]);
                if (params.numInitTrees < 0)
                    throw "-ninit must be non-negative";
				if (params.numInitTrees < params.numNNITrees)
					params.numNNITrees = params.numInitTrees;
				continue;
			}
			if (strcmp(argv[cnt], "-fast") == 0 || strcmp(argv[cnt], "--fast") == 0) {
                // fast search option to resemble FastTree
                if (params.gbo_replicates != 0) {
                    throw("Ultrafast bootstrap (-bb) does not work with -fast option");
                }
                params.numInitTrees = 2;
                if (params.min_iterations == -1)
                    params.min_iterations = 2;
				params.stop_condition = SC_FIXED_ITERATION;
                params.modelEps = 0.05;
                params.suppress_list_of_sequences = true;
                params.suppress_zero_distance_warnings = true;
                params.suppress_duplicate_sequence_warnings = true;
                params.optimize_alg_freerate = "1-BFGS";
                params.opt_gammai = false;
                params.treemix_eps = 0.01;
                params.treemixhmm_eps = 0.05;
                continue;
            }
			if (strcmp(argv[cnt], "-fss") == 0) {
				params.fixStableSplits = true;
//				params.five_plus_five = true;
				continue;
			}
            if (strcmp(argv[cnt], "--stable-thres") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --stable-thres <support_value_threshold>";
                params.stableSplitThreshold = convert_double(argv[cnt]);
                continue;
            }
			if (strcmp(argv[cnt], "-ff") == 0) {
				params.five_plus_five = true;
				continue;
			}

			if (strcmp(argv[cnt], "-tabu") == 0) {
                params.fixStableSplits = true;
				params.tabu = true;
                params.maxCandidates = params.numSupportTrees;
				continue;
			}

            if (strcmp(argv[cnt], "--adt-pert") == 0) {
                if (params.tabu == true) {
                    throw("option -tabu and --adt-pert cannot be combined");
                }
                params.adaptPertubation = true;
                params.stableSplitThreshold = 1.0;
                continue;
            }

            if (strcmp(argv[cnt], "-memcheck") == 0) {
                params.memCheck = true;
                continue;
            }

			if (strcmp(argv[cnt], "--ntop") == 0 || strcmp(argv[cnt], "-ntop") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ntop <number_of_top_parsimony_trees>";
				params.numNNITrees = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "--num-sup-trees") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use --num-sup-trees <number_of_support_trees>";
				params.numSupportTrees = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-fixai") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -fixai <alpha_invar_file>";
				params.alpha_invar_file = argv[cnt];
				continue;
			}

            if (strcmp(argv[cnt], "--opt-gamma-inv") == 0) {
                params.opt_gammai = true;
                continue;
            }

            if (strcmp(argv[cnt], "--no-opt-gamma-inv") == 0) {
                params.opt_gammai = false;
                continue;
            }

            if (strcmp(argv[cnt], "--opt-gammai-fast") == 0) {
                params.opt_gammai_fast = true;
                params.opt_gammai = true;
                continue;
            }

            if (strcmp(argv[cnt], "--opt-gammai-kb") == 0) {
                params.opt_gammai_keep_bran = true;
                params.opt_gammai = true;
                continue;
            }

            if (strcmp(argv[cnt], "--adaptive-eps") == 0) {
                params.testAlphaEpsAdaptive = true;
                continue;
            }
            if (strcmp(argv[cnt], "--rand-alpha") == 0) {
                params.randomAlpha = true;
                continue;
            }

            if (strcmp(argv[cnt], "-eai") == 0) {
                params.exh_ai = true;
                continue;
            }
			if (strcmp(argv[cnt], "-poplim") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -poplim <max_pop_size>";
				params.maxCandidates = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "--nbest") == 0 ||strcmp(argv[cnt], "-nbest") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nbest <number_of_candidate_trees>";
				params.popSize = convert_int(argv[cnt]);
				ASSERT(params.popSize < params.numInitTrees);
				continue;
			}
			if (strcmp(argv[cnt], "-beststart") == 0) {
				params.bestStart = true;
				cnt++;
				if (cnt >= argc)
					throw "Use -best_start <binary_alignment_file>";
				params.binary_aln_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-pll") == 0) {
                throw("-pll option is discontinued.");
				params.pll = true;
				continue;
			}
			if (strcmp(argv[cnt], "-me") == 0 || strcmp(argv[cnt], "--epsilon") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -me <model_epsilon>";
				params.modelEps = convert_double(argv[cnt]);
				if (params.modelEps <= 0.0)
					throw "Model epsilon must be positive";
				if (params.modelEps > 1.0)
					throw "Model epsilon must not be larger than 1.0";
                params.treemix_eps = params.modelEps;
                params.treemixhmm_eps = params.modelEps;
				continue;
			}

            if (strcmp(argv[cnt], "-fundi-epsilon") == 0 || strcmp(argv[cnt], "--fundi-epsilon") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --fundi-epsilon <fundi_epsilon>";
                params.fundiEps = convert_double(argv[cnt]);
                if (params.fundiEps <= 0.0)
                    throw "Fundi epsilon must be positive";
                if (params.fundiEps > 1.0)
                    throw "Fundi epsilon must not be larger than 1.0";
                continue;
            }

            if (strcmp(argv[cnt], "--mf-epsilon") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --mf-epsilon <modelfinder_epsilon>";
                params.modelfinder_eps = convert_double(argv[cnt]);
                if (params.modelfinder_eps <= 0.0)
                    throw "ModelFinder epsilon must be positive";
                if (params.modelfinder_eps > 1.0)
                    throw "ModelFinder epsilon must not be larger than 1.0";
                continue;
            }

            if (strcmp(argv[cnt], "-pars_ins") == 0) {
				params.reinsert_par = true;
				continue;
			}
			if (strcmp(argv[cnt], "-allnni") == 0 || strcmp(argv[cnt], "--allnni") == 0) {
				params.speednni = false;
				continue;
			}
            
			if (strcmp(argv[cnt], "-snni") == 0) {
				params.snni = true;
				// dont need to turn this on here
				//params.autostop = true;
				//params.speednni = true;
				// Minh: why do you turn this on? it doubles curPerStrength at some point
				//params.adaptPert = true;
				continue;
			}
			if (strcmp(argv[cnt], "-iqpnni") == 0) {
				params.snni = false;
				params.start_tree = STT_BIONJ;
				params.numNNITrees = 1;
//            continue; } if (strcmp(argv[cnt], "-auto") == 0) {
//            	params.autostop = true;
				continue;
			}
			if (strcmp(argv[cnt], "--nstop") == 0 || strcmp(argv[cnt], "-nstop") == 0) {
				if (params.stop_condition != SC_BOOTSTRAP_CORRELATION)
					params.stop_condition = SC_UNSUCCESS_ITERATION;
				cnt++;
				if (cnt >= argc)
					throw "Use -nstop <#iterations>";
				params.unsuccess_iteration = convert_int(argv[cnt]);
                if (params.unsuccess_iteration <= 0)
                    throw "-nstop iterations must be positive";
                params.max_iterations = max(params.max_iterations, params.unsuccess_iteration*10);
				continue;
			}
			if (strcmp(argv[cnt], "-lsbran") == 0) {
				params.leastSquareBranch = true;
				continue;
			}
			if (strcmp(argv[cnt], "-manuel") == 0) {
				params.manuel_analytic_approx = true;
				continue;
			}
			if (strcmp(argv[cnt], "-parsbran") == 0) {
				params.pars_branch_length = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bayesbran") == 0) {
				params.bayes_branch_length = true;
				continue;
			}
			if (strcmp(argv[cnt], "-fivebran") == 0
					|| strcmp(argv[cnt], "-nni5") == 0) {
				params.nni5 = true;
				params.nni_type = NNI5;
				continue;
			}
			if (strcmp(argv[cnt], "-onebran") == 0
					|| strcmp(argv[cnt], "-nni1") == 0) {
				params.nni_type = NNI1;
				params.nni5 = false;
				continue;
			}
            
            if (strcmp(argv[cnt], "-nni-eval") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nni-eval <num_evaluation>";
                params.nni5_num_eval = convert_int(argv[cnt]);
                if (params.nni5_num_eval < 1)
                    throw("Positive -nni-eval expected");
                continue;
            }

            if (strcmp(argv[cnt], "-bl-eval") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bl-eval <num_evaluation>";
                params.brlen_num_traversal = convert_int(argv[cnt]);
                if (params.brlen_num_traversal < 1)
                    throw("Positive -bl-eval expected");
                continue;
            }
            
			if (strcmp(argv[cnt], "-smooth") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -smooth <num_iterations>";
				params.numSmoothTree = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-lsnni") == 0) {
				params.leastSquareNNI = true;
				continue;
			}
			if (strcmp(argv[cnt], "-lsvar") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lsvar <o|ft|fm|st|p>";
				if (strcmp(argv[cnt], "o") == 0
						|| strcmp(argv[cnt], "ols") == 0) {
					params.ls_var_type = OLS;
					continue;
				}
				if (strcmp(argv[cnt], "ft") == 0
						|| strcmp(argv[cnt], "first_taylor") == 0) {
					params.ls_var_type = WLS_FIRST_TAYLOR;
					continue;
				}
				if (strcmp(argv[cnt], "fm") == 0
						|| strcmp(argv[cnt], "fitch_margoliash") == 0) {
					params.ls_var_type = WLS_FITCH_MARGOLIASH;
					continue;
				}
				if (strcmp(argv[cnt], "st") == 0
						|| strcmp(argv[cnt], "second_taylor") == 0) {
					params.ls_var_type = WLS_SECOND_TAYLOR;
					continue;
				}
				if (strcmp(argv[cnt], "p") == 0
						|| strcmp(argv[cnt], "pauplin") == 0) {
					params.ls_var_type = WLS_PAUPLIN;
				} else {
					throw "Use -lsvar <o|ft|fm|st|p>";
				}
				continue;
			}
			if (strcmp(argv[cnt], "-eps") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -eps <log-likelihood epsilon>";
				params.loglh_epsilon = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-pb") == 0) { // Enable parsimony branch length estimation
				params.parbran = true;
				continue;
			}
			if (strcmp(argv[cnt], "-x") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -x <iteration_multiple>";
				params.iteration_multiple = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-sp_iter") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sp_iter <number_iteration>";
				params.speedup_iter = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-avh") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -avh <arndt_#bootstrap>";
				params.avh_test = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bootlh") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bootlh <#replicates>";
				params.bootlh_test = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bootpart") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bootpart <part1_length,part2_length,...>";
				params.bootlh_partitions = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-AIC") == 0) {
				params.model_test_criterion = MTC_AIC;
				continue;
			}
			if (strcmp(argv[cnt], "-AICc") == 0 || strcmp(argv[cnt], "-AICC") == 0) {
				params.model_test_criterion = MTC_AICC;
				continue;
			}
			if (strcmp(argv[cnt], "-merit") == 0 || strcmp(argv[cnt], "--merit") == 0) {
                cnt++;
				if (cnt >= argc)
					throw "Use -merit AIC|AICC|BIC";
                if (strcmp(argv[cnt], "AIC") == 0)
                    params.model_test_criterion = MTC_AIC;
                else if (strcmp(argv[cnt], "AICc") == 0 || strcmp(argv[cnt], "AICC") == 0)
                    params.model_test_criterion = MTC_AICC;
                else if (strcmp(argv[cnt], "BIC") == 0)
                    params.model_test_criterion = MTC_BIC;
                else throw "Use -merit AIC|AICC|BIC";
				continue;
			}
			if (strcmp(argv[cnt], "-ms") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ms <model_test_sample_size>";
				params.model_test_sample_size = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-nt") == 0 || strcmp(argv[cnt], "-c") == 0 ||
                strcmp(argv[cnt], "-T") == 0  || strcmp(argv[cnt], "--threads") == 0) {
				cnt++;
				if (cnt >= argc)
				throw "Use -nt <num_threads|AUTO>";
                if (iEquals(argv[cnt], "AUTO"))
                    params.num_threads = 0;
                else {
                    params.num_threads = convert_int(argv[cnt]);
                    if (params.num_threads < 1)
                        throw "At least 1 thread please";
                }
				continue;
			}
            
            if (strcmp(argv[cnt], "-ntmax") == 0 || strcmp(argv[cnt], "--threads-max") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -ntmax <num_threads_max>";
                params.num_threads_max = convert_int(argv[cnt]);
                if (params.num_threads_max < 1)
                    throw "At least 1 thread please";
                continue;
            }
            
            if (strcmp(argv[cnt], "--thread-model") == 0) {
                params.openmp_by_model = true;
                continue;
            }

            if (strcmp(argv[cnt], "--thread-site") == 0) {
                params.openmp_by_model = false;
                continue;
            }

//			if (strcmp(argv[cnt], "-rootstate") == 0) {
//                cnt++;
//                if (cnt >= argc)
//                    throw "Use -rootstate <rootstate>";
//                params.root_state = argv[cnt];
//                params.SSE = LK_NORMAL;
//                continue;
//			}
			if (strcmp(argv[cnt], "-ct") == 0) {
            	params.count_trees = true;
            	continue;
			}
			if (strcmp(argv[cnt], "--sprrad") == 0 || strcmp(argv[cnt], "--radius") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sprrad <SPR radius used in parsimony search>";
				params.sprDist = convert_int(argv[cnt]);
				continue;
			}
            
            if (strcmp(argv[cnt], "--mpcost") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --mpcost <parsimony_cost_file>";
                params.sankoff_cost_file = argv[cnt];
                continue;
            }
            
			if (strcmp(argv[cnt], "-no_rescale_gamma_invar") == 0) {
				params.no_rescale_gamma_invar = true;
				continue;
			}

			if (strcmp(argv[cnt], "-wsi") == 0) {
				params.compute_seq_identity_along_tree = true;
				continue;
			}
            
            if (strcmp(argv[cnt], "--no-seq-comp") == 0) {
                params.compute_seq_composition = false;
                continue;
            }
            
			if (strcmp(argv[cnt], "-t") == 0 || strcmp(argv[cnt], "-te") == 0 || strcmp(argv[cnt], "--tree") == 0) {
                if (strcmp(argv[cnt], "-te") == 0) {
                    if (params.gbo_replicates != 0) {
                        throw("Ultrafast bootstrap does not work with -te option");
                    }
                    params.min_iterations = 0;
                    params.stop_condition = SC_FIXED_ITERATION;
                }
				cnt++;
				if (cnt >= argc)
					throw "Use -t,-te <start_tree | BIONJ | PARS | PLLPARS | RANDOM>";
				else if (strcmp(argv[cnt], "PARS") == 0)
					params.start_tree = STT_PARSIMONY;
				else if (strcmp(argv[cnt], "PLLPARS") == 0)
					params.start_tree = STT_PLL_PARSIMONY;
                else if (strcmp(argv[cnt], "RANDOM") == 0 || strcmp(argv[cnt], "RAND") == 0)
                {
					params.start_tree = STT_RANDOM_TREE;
                    params.tree_gen = YULE_HARDING;
                    
                    // <NUM_TAXA> is required if an alignment is not provided
                    if ((params.original_params.find("-s ") == std::string::npos)
                        && (params.original_params.find("--aln ") == std::string::npos))
                        throw "Use -t RANDOM{<MODEL>,<NUM_TAXA>} where <MODEL> is yh, u, cat, bal, bd{<birth_rate>,<death_rate>} stands for Yule-Harding, Uniform, Caterpillar, Balanced, Birth-Death model respectively; <NUM_TAXA> could be a fixed number, a list (<NUM_1>,<NUM_2>,...,<NUM_N>), or a Uniform distribution U(<LOWER_BOUND>,<UPPER_BOUND>). <NUM_TAXA> is only required if an alignment is not provided by -s <ALIGNMENT_FILE>.";
                }
                else if (START_TREE_RECOGNIZED(argv[cnt])) {
                    params.start_tree_subtype_name = argv[cnt];
                    params.start_tree = STT_BIONJ;
                }
                else
                {
                    string ERR_MSG = "Use -t RANDOM{<MODEL>,<NUM_TAXA>} where <MODEL> is yh, u, cat, bal, bd{<birth_rate>,<death_rate>} stands for Yule-Harding, Uniform, Caterpillar, Balanced, Birth-Death model respectively; <NUM_TAXA> could be a fixed number, a list (<NUM_1>,<NUM_2>,...,<NUM_N>), or a Uniform distribution U(<LOWER_BOUND>,<UPPER_BOUND>). <NUM_TAXA> is only required if an alignment is not provided by -s <ALIGNMENT_FILE>.";
                    string t_params = argv[cnt];
                    transform(t_params.begin(), t_params.end(), t_params.begin(), ::toupper);
                    string KEYWORD = "RANDOM";
                    if ((t_params.length() > KEYWORD.length())
                        && (!t_params.substr(0, KEYWORD.length()).compare(KEYWORD)))
                    {
                        params.start_tree = STT_RANDOM_TREE;
                        
                        // validate the input
                        if ((t_params[KEYWORD.length()]!='{')
                            ||(t_params[t_params.length()-1]!='}'))
                        {
                            // try to parse the number of taxa -> users may input RANDOM<NUM> here
                            string num_taxa_str = t_params.substr(KEYWORD.length(), t_params.length() - KEYWORD.length());
                            // check if num_taxa_str presents an integer
                            if (std::all_of(num_taxa_str.begin(), num_taxa_str.end(), ::isdigit))
                            {
                                // set the default model is yule harding
                                params.tree_gen = YULE_HARDING;

                                // set the number of taxa
                                params.sub_size = convert_int(num_taxa_str.c_str());

                                continue;
                            }
                            else
                                throw ERR_MSG;
                        }
                        
                        // remove "RANDOM{"
                        t_params.erase(0, KEYWORD.length() + 1);
                        
                        // remove "}"
                        t_params = t_params.substr(0, t_params.length()-1);
                        
                        // extract model_name & #taxa
                        // detect delimiter
                        string delimiter = ",";
                        bool inside_bracket = false;
                        for (int i = 0; i < t_params.length(); i++)
                        {
                            // detect brackets
                            if (t_params[i] == '{')
                                inside_bracket = true;
                            else if (t_params[i] == '}')
                                inside_bracket = false;
                               
                            // only check separator outside brackets
                            if (!inside_bracket && (t_params[i] == ',' || t_params[i] == '/'))
                            {
                                delimiter = t_params[i];
                                break;
                            }
                        }
                        
                        // change the delimiter to "}," in case with birth-death model
                        if ((t_params.length() > 2)
                            && (!t_params.substr(0, 2).compare("BD")))
                            delimiter = "}" + delimiter;
                        
                        string model_name = t_params;
                        string num_taxa = "";
                        if (t_params.find(delimiter) != string::npos)
                        {
                            int delimiter_index = t_params.find(delimiter);
                            model_name = t_params.substr(0, delimiter_index +  delimiter.length() - 1);
                            num_taxa = t_params;
                            num_taxa.erase(0, delimiter_index + delimiter.length());
                        }
                        
                        // <NUM_TAXA> is required if an alignment is not provided
                        if ((params.original_params.find("-s ") == std::string::npos)
                            && (params.original_params.find("--aln ") == std::string::npos)
                            && num_taxa.length() == 0)
                            throw ERR_MSG;
                        
                        // parse #taxa
                        if (num_taxa.length() > 0)
                        {
                            // handle the case with a list (<NUM_1>,<NUM_2>,...,<NUM_N>)
                            if (num_taxa[0] == '{')
                            {
                                // detect the seperator
                                delimiter = ",";
                                if (num_taxa.find('/') != std::string::npos)
                                    delimiter = "/";
                                
                                // validate the input
                                if ((num_taxa[0]!='{')
                                    ||(num_taxa[num_taxa.length()-1]!='}'))
                                    throw ERR_MSG;
                                
                                // remove "U{"
                                num_taxa.erase(0, 1);
                                
                                // remove "}"
                                num_taxa = num_taxa.substr(0, num_taxa.length()-1);
                                
                                // get the list of numbers of taxa
                                while (num_taxa.length()>0)
                                {
                                    // get the <NUM_X>
                                    int tmp_num = convert_int(num_taxa.substr(0, num_taxa.find(delimiter)).c_str());
                                    if (tmp_num <= 3)
                                        throw "<NUM_TAXA> must be greater than 3.";
                                    params.alisim_num_taxa_list.push_back(tmp_num);
                                    
                                    // remove "<NUM_TAXA>,"
                                    if (num_taxa.find(delimiter) != string::npos)
                                        num_taxa.erase(0, num_taxa.find(delimiter) + delimiter.length());
                                    else
                                        num_taxa = "";
                                }
                            }
                            // handle the case with a Uniform Distribution U(<LOWER_BOUND>,<UPPER_BOUND>)
                            else if (num_taxa[0] == 'U')
                            {
                                // detect the seperator
                                delimiter = ",";
                                if (num_taxa.find('/') != std::string::npos)
                                    delimiter = "/";
                                
                                // validate the input
                                if ((num_taxa[1]!='{')
                                    ||(num_taxa[num_taxa.length()-1]!='}')
                                    || (num_taxa.find(delimiter) == string::npos))
                                    throw ERR_MSG;
                                
                                // remove "U{"
                                num_taxa.erase(0, 2);
                                
                                // remove "}"
                                num_taxa = num_taxa.substr(0, num_taxa.length()-1);
                                
                                // get the <LOWER_BOUND>
                                params.alisim_num_taxa_uniform_start = convert_int(num_taxa.substr(0, num_taxa.find(delimiter)).c_str());
                                if (params.alisim_num_taxa_uniform_start <= 3)
                                    throw "<LOWER_BOUND> must be greater than 3.";
                                
                                // remove "<LOWER_BOUND>,"
                                num_taxa.erase(0, num_taxa.find(delimiter) + delimiter.length());
                                
                                // get <UPPER_BOUND>
                                params.alisim_num_taxa_uniform_end = convert_int(num_taxa.c_str());
                                if (params.alisim_num_taxa_uniform_start > params.alisim_num_taxa_uniform_end)
                                    throw "<UPPER_BOUND> must not be less than <LOWER_BOUND>";
                                
                            }
                            // handle the case with a fixed number
                            else
                            {
                                params.sub_size = convert_int(num_taxa.c_str());
                                if (params.sub_size <= 3)
                                    throw ERR_MSG +" <NUM_TAXA> must be greater than 3.";
                            }
                        }
                        
                        transform(model_name.begin(), model_name.end(), model_name.begin(), ::toupper);
                        // parse model
                        if (model_name == "YH")
                            params.tree_gen = YULE_HARDING;
                        else if (model_name == "U")
                            params.tree_gen = UNIFORM;
                        else if (model_name == "CAT")
                            params.tree_gen = CATERPILLAR;
                        else if (model_name == "BAL")
                            params.tree_gen = BALANCED;
                        else
                        {
                            KEYWORD = "BD";
                            if ((model_name.length() > KEYWORD.length())
                                && (!model_name.substr(0, KEYWORD.length()).compare(KEYWORD)))
                            {
                                params.tree_gen = BIRTH_DEATH;
                                
                                // detect the seperator
                                delimiter = ",";
                                if (model_name.find('/') != std::string::npos)
                                    delimiter = "/";
                                
                                // validate the input
                                if ((model_name[KEYWORD.length()]!='{')
                                    ||(model_name[model_name.length()-1]!='}')
                                    ||(model_name.find(delimiter) == string::npos))
                                    throw "Use bd{<birth_rate>,<death_rate>} to specify the Birth-Death model.";
                                
                                // remove "bd{"
                                model_name.erase(0, KEYWORD.length() + 1);
                                
                                // remove "}"
                                model_name = model_name.substr(0, model_name.length()-1);
                                
                                // get birth_rate
                                params.birth_rate = convert_double(model_name.substr(0, model_name.find(delimiter)).c_str());
                                if (params.birth_rate <= 0)
                                    throw "<birth_rate> must be positive.";
                                
                                // remove "<birth_rate>,"
                                model_name.erase(0, model_name.find(delimiter) + delimiter.length());
                                
                                // get death_rate
                                params.death_rate = convert_double(model_name.c_str());
                                if (params.death_rate < 0 || params.death_rate >= params.birth_rate)
                                    throw "<death_rate> must be non-negative and less than <birth_rate>";
                                
                                // normalize birth_rate & death_rate
                                double sum_rate = params.death_rate + params.birth_rate;
                                if (fabs(sum_rate-1.0) > 1e-5)
                                {
                                    outWarning("Normalizing birth rate and death rate so that sum of them is equal to 1.");
                                    params.death_rate /= sum_rate;
                                    params.birth_rate /= sum_rate;
                                }
                            }
                            else
                                throw ERR_MSG;
                        }
                    }
                    else {
                        params.user_file = argv[cnt];
                        if (params.min_iterations == 0)
                            params.start_tree = STT_USER_TREE;
                    }
                }
				continue;
			}
            
            if (strcmp(argv[cnt], "--no-ml-tree") == 0) {
                params.modelfinder_ml_tree = false;
                continue;
            }
            
            if (strcmp(argv[cnt], "--tree-fix") == 0) {
                if (params.gbo_replicates != 0) {
                    outError("Ultrafast bootstrap does not work with -te option");
                }
                params.min_iterations = 0;
                params.stop_condition = SC_FIXED_ITERATION;
                continue;
            }
                
            if (strcmp(argv[cnt], "-g") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -g <constraint_tree>";
                params.constraint_tree_file = argv[cnt];
                params.ignore_identical_seqs = false;
                continue;
            }
            
			if (strcmp(argv[cnt], "-lmap") == 0 || strcmp(argv[cnt], "--lmap") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lmap <likelihood_mapping_num_quartets>";
                if (iEquals(argv[cnt], "ALL")) {
                    params.lmap_num_quartets = 0;
                } else {
                    params.lmap_num_quartets = convert_int64(argv[cnt]);
                    if (params.lmap_num_quartets < 0)
                        throw "Number of quartets must be >= 1";
                }
				continue;
			}

			if (strcmp(argv[cnt], "-lmclust") == 0 || strcmp(argv[cnt], "--lmclust") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lmclust <likelihood_mapping_cluster_file>";
				params.lmap_cluster_file = argv[cnt];
				// '-keep_ident' is currently required to allow a 1-to-1 mapping of the 
				// user-given groups (HAS) - possibly obsolete in the future versions
				params.ignore_identical_seqs = false;
                if (params.lmap_num_quartets < 0)
                    params.lmap_num_quartets = 0;
				continue;
			}

			if (strcmp(argv[cnt], "-wql") == 0 || strcmp(argv[cnt], "--quartetlh") == 0) {
				params.print_lmap_quartet_lh = true;
				continue;
			}

			if (strcmp(argv[cnt], "-mixlen") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mixlen <number of mixture branch lengths for heterotachy model>";
				params.num_mixlen = convert_int(argv[cnt]);
				if (params.num_mixlen < 1)
					throw("-mixlen must be >= 1");
				continue;
			}
            
			if (strcmp(argv[cnt], "--link-alpha") == 0) {
				params.link_alpha = true;
				continue;
			}

            if (strcmp(argv[cnt], "--link-model") == 0) {
                params.link_model = true;
                continue;
            }

            if (strcmp(argv[cnt], "--model-joint") == 0 || strcmp(argv[cnt], "--link-partition") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use " + string(argv[cnt-1]) + " MODEL_NAME";
                params.model_joint = argv[cnt];
                params.link_model = true;
                continue;
            }

            if (strcmp(argv[cnt], "--unlink-tree") == 0) {
                params.partition_type = TOPO_UNLINKED;
                params.ignore_identical_seqs = false;
                continue;
            }
            
			if (strcmp(argv[cnt], "-redo") == 0 || strcmp(argv[cnt], "--redo") == 0) {
				params.ignore_checkpoint = true;
                // 2020-04-27: SEMANTIC CHANGE: also redo ModelFinder
                params.model_test_again = true;
				continue;
			}

            if (strcmp(argv[cnt], "-tredo") == 0 || strcmp(argv[cnt], "--tredo") == 0 || strcmp(argv[cnt], "--redo-tree") == 0) {
                params.ignore_checkpoint = true;
                continue;
            }

			if (strcmp(argv[cnt], "-undo") == 0 || strcmp(argv[cnt], "--undo") == 0) {
				params.force_unfinished = true;
				continue;
			}

			if (strcmp(argv[cnt], "-cptime") == 0 || strcmp(argv[cnt], "--cptime") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -cptime <checkpoint_time_interval>";
				params.checkpoint_dump_interval = convert_int(argv[cnt]);
				continue;
			}
            
            if (strcmp(argv[cnt], "--all-checkpoint") == 0) {
                params.print_all_checkpoints = true;
                continue;
            }
            
			if (strcmp(argv[cnt], "--no-log") == 0) {
				params.suppress_output_flags |= OUT_LOG;
				continue;
			}

			if (strcmp(argv[cnt], "--no-treefile") == 0) {
				params.suppress_output_flags |= OUT_TREEFILE;
				continue;
			}
			if (strcmp(argv[cnt], "--no-iqtree") == 0) {
				params.suppress_output_flags |= OUT_IQTREE;
				continue;
			}
			if (strcmp(argv[cnt], "--no-outfiles") == 0) {
				params.suppress_output_flags |= OUT_LOG + OUT_TREEFILE + OUT_IQTREE;
				continue;
			}

            // -- Mon Apr 17 21:18:23 BST 2017
            // DONE Minh: merged correctly.
            if (strcmp(argv[cnt], "--scaling-squaring") == 0) {
                params.matrix_exp_technique = MET_SCALING_SQUARING;
                continue;
            }
            if (strcmp(argv[cnt], "--eigenlib") == 0) {
                params.matrix_exp_technique = MET_EIGEN3LIB_DECOMPOSITION;
                continue;
            }
            if (strcmp(argv[cnt], "--eigen") == 0) {
                params.matrix_exp_technique = MET_EIGEN_DECOMPOSITION;
                continue;
            }
            if (strcmp(argv[cnt], "--lie-markov") == 0) {
                params.matrix_exp_technique = MET_LIE_MARKOV_DECOMPOSITION;
                continue;
            }            
			if (strcmp(argv[cnt], "--no-uniqueseq") == 0) {
				params.suppress_output_flags |= OUT_UNIQUESEQ;
				continue;
			}
            // --

            if (strcmp(argv[cnt], "--dating") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --dating LSD or --dating mcmctree";
                params.dating_method = argv[cnt];
                if (params.dating_method != "LSD"){
                    if (params.dating_method != "mcmctree"){
                        throw "Use --dating LSD or --dating mcmctree";
                    }

                }
                continue;
            }

            if (strcmp(argv[cnt], "--date") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --date <date_file>|TAXNAME";
                if (params.dating_method == "")
                    params.dating_method = "LSD";
                params.date_file = argv[cnt];
                continue;
            }

            if (strcmp(argv[cnt], "--date-tip") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --date-tip <YYYY[-MM-DD]>";
                if (params.dating_method == "")
                    params.dating_method = "LSD";
                params.date_tip = argv[cnt];
                continue;
            }

            if (strcmp(argv[cnt], "--date-root") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --date-root <YYYY[-MM-DD]>";
                if (params.dating_method == "")
                    params.dating_method = "LSD";
                params.date_root = argv[cnt];
                continue;
            }

            if (strcmp(argv[cnt], "--date-no-outgroup") == 0) {
                params.date_with_outgroup = false;
                continue;
            }

            if (strcmp(argv[cnt], "--date-outgroup") == 0) {
                params.date_with_outgroup = true;
                continue;
            }

            if (strcmp(argv[cnt], "--date-ci") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --date-ci <number_of_replicates>";
                params.date_replicates = convert_int(argv[cnt]);
                if (params.date_replicates < 1)
                    throw "--date-ci must be positive";
                continue;
            }

            if (strcmp(argv[cnt], "--clock-sd") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --clock-sd <standard_dev_of_lognormal_relaxed_lock>";
                params.clock_stddev = convert_double(argv[cnt]);
                if (params.clock_stddev < 0)
                    throw "--clock-sd must be non-negative";
                continue;
            }

            if (strcmp(argv[cnt], "--date-outlier") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --date-outlier <z_score_for_removing_outlier_nodes>";
                params.date_outlier = convert_double(argv[cnt]);
                if (params.date_outlier < 0)
                    throw "--date-outlier must be non-negative";
                continue;
            }

            if (strcmp(argv[cnt], "--date-debug") == 0) {
                params.date_debug = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "--suppress-list-of-sequences") == 0) {
                params.suppress_list_of_sequences = true;
                continue;
            }

            if (strcmp(argv[cnt], "--suppress-zero-distance") == 0) {
                params.suppress_zero_distance_warnings = true;
                continue;
            }

            if (strcmp(argv[cnt], "--suppress-duplicate-sequence") == 0) {
                params.suppress_duplicate_sequence_warnings = true;
                continue;
            }

            if (strcmp(argv[cnt], "--date-options") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --date-options <extra_options_for_dating_method>";
                params.dating_options = argv[cnt];
                continue;
            }

            if (strcmp(argv[cnt], "--mcmc-clock") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --mcmc-clock <EQUAL|IND|CORR>";
                if (strcmp(argv[cnt], "EQUAL")==0)
                {
                    params.mcmc_clock = EQUAL_RATES;
                }
                else if (strcmp(argv[cnt], "IND")==0)
                {
                    params.mcmc_clock = INDEPENDENT;
                }
                else if (strcmp(argv[cnt], "CORR")==0)
                {
                    params.mcmc_clock = CORRELATED;
                }else
                {
                    throw "Only equal rate, independent and correlated clock models are supported in MCMCtree";
                }
                continue;
            }

            if (strcmp(argv[cnt], "--mcmc-bds") == 0) {
                cnt++;
                params.mcmc_bds = argv[cnt];
                StrVector mcmc_bds_vec;
                convert_string_vec(params.mcmc_bds.c_str(), mcmc_bds_vec, ',');
                if (mcmc_bds_vec.size() != 3 || mcmc_bds_vec[0].empty() ||
                    mcmc_bds_vec[1].empty() ||
                    mcmc_bds_vec[2].empty())
                {
                    throw
                        "three parameters should be set for birth-death model of MCMCtree (birth-rate, death-rate and sampling-fraction)";
                }
                params.mcmc_bds = mcmc_bds_vec[0] + " " + mcmc_bds_vec[1] + " " + mcmc_bds_vec[2];
                continue;
            }

            if (strcmp(argv[cnt], "--mcmc-iter") == 0) {
                cnt++;
                params.mcmc_iter = argv[cnt];
                StrVector mcmc_iter_vec;
                convert_string_vec(params.mcmc_iter.c_str(), mcmc_iter_vec, ',');
                if (mcmc_iter_vec.size() != 3 || mcmc_iter_vec[0].empty() ||
                    mcmc_iter_vec[1].empty() ||
                    mcmc_iter_vec[2].empty())
                {
                    throw "three parameters should be set for MCMCtree dating (Burin, samplefreq and nsamples)";
                }
                continue;
            }

            // added by TD
            // todo: give it a fancy name
            if (strcmp(argv[cnt], "--use-nn-model") == 0) {
                params.use_nn_model = true;
                continue;
            }

            // added by TD
            if (strcmp(argv[cnt], "--nn-path-model") == 0) {
                cnt++;
                params.nn_path_model = argv[cnt];
                params.use_nn_model = true;
                continue;
            }

            // added by TD
            if (strcmp(argv[cnt], "--nn-path-rates") == 0) {
                cnt++;
                params.nn_path_rates = argv[cnt];
                params.use_nn_model = true;
                continue;
            }

            if (arg=="-progress-bar" || arg=="--progress-bar" || arg=="-bar") {
                progress_display::setProgressDisplay(true);
                continue;
            }
            
            if (strcmp(argv[cnt], "--alisim") == 0) {
                params.alisim_active = true;
                params.multi_rstreams_used = true;
                
                cnt++;
                if (cnt >= argc || argv[cnt][0] == '-')
                {
                    usage_alisim();
                    throw "";
                }
                params.alisim_output_filename = argv[cnt];
                
                // set default model_name
                if (params.model_name.empty())
                    params.model_name = "HKY";
                
                continue;
            }

            if (strcmp(argv[cnt], "--index-from-one") == 0) {
                params.site_starting_index = 1;

                continue;
            }
            
            if (strcmp(argv[cnt], "--distribution") == 0) {
                cnt++;
                if (cnt >= argc || argv[cnt][0] == '-')
                    throw "Use --distribution <distribution_file>";
                
                params.alisim_distribution_definitions = argv[cnt];

                continue;
            }
            
            if (strcmp(argv[cnt], "--no-copy-gaps") == 0) {
                params.alisim_no_copy_gaps = true;
                
                continue;
            }
            
            if (strcmp(argv[cnt], "--single-output") == 0) {
                params.alisim_single_output = true;
                
                continue;
            }
            
            if (strcmp(argv[cnt], "--length") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --length <SEQUENCE_LENGTH>";
                params.alisim_sequence_length = convert_int(argv[cnt]);
                if (params.alisim_sequence_length < 1)
                    throw "Positive --length please";
                continue;
            }
            
            if (strcmp(argv[cnt], "--rebuild-indel-history") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --rebuild-indel-history <proportion>";
                params.rebuild_indel_history_param = convert_double(argv[cnt]);
                if (params.rebuild_indel_history_param < 0 || params.rebuild_indel_history_param > 1)
                    throw "<proportion> must be between 0 and 1.";
                continue;
            }
            
            if (strcmp(argv[cnt], "--keep-seq-order") == 0) {
                params.keep_seq_order = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "--mem-limit") == 0) {
                cnt++;
                if (cnt >= argc || argv[cnt][0] == '-')
                    throw "Use --mem-limit <FACTOR>";
                
                params.mem_limit_factor = convert_double(argv[cnt]);
                if (params.mem_limit_factor <= 0 || params.mem_limit_factor > 1)
                    throw "<FACTOR> must be in range (0,1]";

                continue;
            }
            
            if (strcmp(argv[cnt], "--delete-output") == 0) {
                params.delete_output = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "--openmp-alg") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --openmp-alg <ALG>";
                
                string openmp_algorithm = argv[cnt];
                transform(openmp_algorithm.begin(), openmp_algorithm.end(), openmp_algorithm.begin(), ::toupper);
                
                if (openmp_algorithm == "IM")
                    params.alisim_openmp_alg = IM;
                else if (openmp_algorithm == "EM")
                    params.alisim_openmp_alg = EM;
                else
                    throw "AliSim-OpenMP algorithm should be IM or EM.";
                continue;
            }
            
            if (strcmp(argv[cnt], "--no-merge") == 0) {
                params.no_merge = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "--indel-rate-variation") == 0) {
                params.indel_rate_variation = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "--num-alignments") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --num-alignments <NUMBER_OF_DATASETS>";
                params.alisim_dataset_num = convert_int(argv[cnt]);
                if (params.alisim_dataset_num < 1)
                    throw "Positive --num-alignments please";
                continue;
            }
            
            if (strcmp(argv[cnt], "--root-seq") == 0  || strcmp(argv[cnt], "--ref-seq") == 0 || strcmp(argv[cnt], "-ref") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --root-seq <ALN_FILE>,<SEQ_NAME>";
                string ancestral_sequence_params = argv[cnt];
                string delimiter = ",";
                if ((ancestral_sequence_params.find(delimiter) == std::string::npos)||ancestral_sequence_params.find(delimiter) == 0)
                    throw "Use --root-seq <ALN_FILE>,<SEQ_NAME>";
                auto delimiter_pos = ancestral_sequence_params.find(delimiter);
                params.root_ref_seq_aln = ancestral_sequence_params.substr(0, delimiter_pos);
                ancestral_sequence_params.erase(0, delimiter_pos + delimiter.length());
                params.root_ref_seq_name = ancestral_sequence_params;
                if (!params.root_ref_seq_aln.length() || !params.root_ref_seq_name.length())
                    throw "Use --root-seq <ALN_FILE>,<SEQ_NAME>";
                
                continue;
            }
            if (strcmp(argv[cnt], "-aln-format") == 0 || strcmp(argv[cnt], "--aln-format") == 0) {
                cnt++;
                if (cnt >= argc)
                    outError("Use -aln-format MAPLE, PHYLIP, FASTA, or AUTO");
                params.in_aln_format_str = argv[cnt];

                continue;
            }
            if (strcmp(argv[cnt], "--tree-search") == 0 || strcmp(argv[cnt], "-tree-search") == 0) {
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -tree-search <FAST|NORMAL|EXHAUSTIVE>");

                params.tree_search_type_str = argv[cnt];
                continue;
            }
            if (strcmp(argv[cnt], "--shallow-tree-search") == 0 || strcmp(argv[cnt], "-shallow-search") == 0) {

                params.shallow_tree_search = true;

            if (strcmp(argv[cnt], "--mutation") == 0 || strcmp(argv[cnt], "-mut") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --mutation <MUTATION_FILE>";
                params.mutation_file = argv[cnt];
                params.include_pre_mutations = true;
                continue;
            }

                continue;
            }
            if (strcmp(argv[cnt], "--replace-input-tree") == 0 || strcmp(argv[cnt], "-rep-tree") == 0) {
                params.allow_replace_input_tree = true;
                continue;
            }
            if (strcmp(argv[cnt], "--output-multifurcating-tree") == 0 || strcmp(argv[cnt], "-out-mul-tree") == 0) {

                params.tree_format_str = "MUL";

                continue;
            }
            if (strcmp(argv[cnt], "--make-consistent") == 0 || strcmp(argv[cnt], "-consistent") == 0) {
                params.make_consistent = true;
                continue;
            }

            if (argv[cnt][0] == '-') {
                string err = "Invalid \"";
                err += argv[cnt];
                err += "\" option.";
                throw err;
            } else {
                if (params.user_file == nullptr)
                    params.user_file = argv[cnt];
                else
                    params.out_file = argv[cnt];
            }
        }
        // try
        catch (const char *str) {
            if (MPIHelper::getInstance().isMaster())
                outError(str);
            else
                exit(EXIT_SUCCESS);
            //} catch (char *str) {
            //outError(str);
        } catch (string str) {
            if (MPIHelper::getInstance().isMaster())
                outError(str);
            else
                exit(EXIT_SUCCESS);
        } catch (...) {
            string err = "Unknown argument \"";
            err += argv[cnt];
            err += "\"";
            if (MPIHelper::getInstance().isMaster())
                outError(err);
            else
                exit(EXIT_SUCCESS);
        }

    } // for
    if (!params.user_file && !params.aln_file && !params.ngs_file && !params.ngs_mapped_reads && !params.partition_file && !params.alisim_active) {
#ifdef IQ_TREE
        quickStartGuide();
//        usage_iqtree(argv, false);
#else
        usage(argv, false);
#endif
    }

    if (params.treemix_optimize_methods.find("hmm")!=string::npos &&
        params.model_name.find("+T") != string::npos) {
        params.optimize_params_use_hmm = true;
    } else {
        params.optimize_params_use_hmm = false;
    }

//    if (params.do_au_test)
//        outError("The AU test is temporarily disabled due to numerical issue when bp-RELL=0");

    if (params.root != nullptr && params.is_rooted)
        outError("Not allowed to specify both -o <taxon> and -root");
    
    if (params.model_test_and_tree && params.partition_type != BRLEN_OPTIMIZE)
        outError("-mtree not allowed with edge-linked partition model (-spp or -q)");
    
    if (params.do_au_test && params.topotest_replicates == 0)
        outError("For AU test please specify number of bootstrap replicates via -zb option");
    
    if (params.lh_mem_save == LM_MEM_SAVE && params.partition_file)
        outError("-mem option does not work with partition models yet");
    
    if (params.gbo_replicates && params.num_bootstrap_samples)
        outError("UFBoot (-bb) and standard bootstrap (-b) must not be specified together");
    
    if ((params.model_name.find("ONLY") != string::npos || (params.model_name.substr(0,2) == "MF" && params.model_name.substr(0,3) != "MFP")) && (params.gbo_replicates || params.num_bootstrap_samples))
        outError("ModelFinder only cannot be combined with bootstrap analysis");
    
    if (params.num_runs > 1 && !params.treeset_file.empty())
        outError("Can't combine --runs and -z options");
    
    if (params.num_runs > 1 && params.lmap_num_quartets >= 0)
        outError("Can't combine --runs and -lmap options");

    if (params.terrace_analysis_tphast && !params.partition_file)
        params.terrace_analysis_tphast = false;

    if (params.constraint_tree_file && params.partition_type == TOPO_UNLINKED)
        outError("-g constraint tree option does not work with -S yet.");

    if (params.num_bootstrap_samples && params.partition_type == TOPO_UNLINKED)
        outError("-b bootstrap option does not work with -S yet.");

    //added to remove situations where we're optimizing a linked rate matrix when we really shouldn't be -JD
    if (params.optimize_linked_gtr && params.model_name.find("GTR") == string::npos && params.model_joint.find("GTR") == string::npos)
        outError("Must have either GTR or GTR20 as part of the model when using --link-exchange-rates.");

    if (params.use_nn_model && params.modelomatic)
        outError("--modelomatic option does not work with --use-nn-model.");

    if (params.dating_method != "") {
    #ifndef USE_LSD2
        outError("IQ-TREE was not compiled with LSD2 library, rerun cmake with -DUSE_LSD2=ON option");
    #endif
    }

    if (params.date_file.empty()) {
        if (params.date_root.empty() ^ params.date_tip.empty())
            outError("Both --date-root and --date-tip must be provided when --date file is absent");
    }
    
	// Diep:
	if(params.ufboot2corr == true){
		if(params.gbo_replicates <= 0) params.ufboot2corr = false;
		else params.stop_condition = SC_UNSUCCESS_ITERATION;

		params.print_ufboot_trees = 2; // 2017-09-25: fix bug regarding the order of -bb 1000 -bnni -wbt
	}

    if (!params.out_prefix) {
    	if (params.eco_dag_file)
    		params.out_prefix = params.eco_dag_file;
        else if (params.user_file && params.consensus_type == CT_ASSIGN_SUPPORT_EXTENDED)
            params.out_prefix = params.user_file;
        else if (params.partition_file) {
            params.out_prefix = params.partition_file;
            if (params.out_prefix[strlen(params.out_prefix)-1] == '/' || params.out_prefix[strlen(params.out_prefix)-1] == '\\') {
                params.out_prefix[strlen(params.out_prefix)-1] = 0;
            }
        } else if (params.aln_file) {
            params.out_prefix = params.aln_file;
            if (params.out_prefix[strlen(params.out_prefix)-1] == '/' || params.out_prefix[strlen(params.out_prefix)-1] == '\\') {
                params.out_prefix[strlen(params.out_prefix)-1] = 0;
            }
        } else if (params.ngs_file)
            params.out_prefix = params.ngs_file;
        else if (params.ngs_mapped_reads)
            params.out_prefix = params.ngs_mapped_reads;
        else
            params.out_prefix = params.user_file;
    }

    if (params.model_name.find("LINK") != string::npos || params.model_name.find("MERGE") != string::npos)
        if (params.partition_merge == MERGE_NONE)
            params.partition_merge = MERGE_RCLUSTERF;
    
    if (params.alisim_active && !params.aln_file && !params.user_file && !params.partition_file && params.tree_gen == NONE)
        outError("A tree filepath is a mandatory input to execute AliSim when neither Inference mode nor Random mode (generating a random tree) is inactive. Use -t <TREE_FILEPATH> ; or Activate the inference mode by -s <ALIGNMENT_FILE> ; or Activate Random mode by -t RANDOM{<MODEL>,<NUM_TAXA>} where <MODEL> is yh, u, cat, bal, bd{<birth_rate>,<death_rate>} stands for Yule-Harding, Uniform, Caterpillar, Balanced, Birth-Death model respectively.");
    // terminate if using AliSim with -ft or -fs site-specific model (ModelSet)
    // computeTransMatix has not yet implemented for ModelSet
    if (params.alisim_active && (params.tree_freq_file || params.site_freq_file))
        outError("Sorry! `-ft` (--site-freq) and `-fs` (--tree-freq) options are not fully supported in AliSim. However, AliSim can estimate posterior mean frequencies from the alignment. Please try again without `-ft` and `-fs` options!");
    
    // set default filename for the random tree if AliSim is running in Random mode
    if (params.alisim_active && !params.user_file && params.tree_gen != NONE)
    {
        if (!params.aln_file)
        {
            string tmp_user_file = params.alisim_output_filename+".treefile";
            params.user_file = new char[tmp_user_file.length() + 1];
            strcpy(params.user_file, tmp_user_file.c_str());
        }
        else
        {
            // initialize output_filepath
            std::string output_filepath(params.aln_file);
            output_filepath = output_filepath.substr(0, output_filepath.find_last_of("/\\") + 1);
            output_filepath = output_filepath + params.alisim_output_filename +".treefile";
            
            params.user_file = new char[output_filepath.length() + 1];
            strcpy(params.user_file, output_filepath.c_str());
        }
        params.out_prefix = params.user_file;
    }
    
    // set the number of sites per state at 3 (for CODON)
    if (params.sequence_type)
    {
        std::string sequence_type(params.sequence_type);
    }

    // parse the profile mixture model
    // R+Fx -> MIX{S+FO,S+FO,...,S+FO} with x classes and S is a linked substitution matrix (i.e. linked exchangeabilities)
    parseProfileMix(params);

    //    if (MPIHelper::getInstance().isWorker()) {
    // BUG: setting out_prefix this way cause access to stack, which is cleaned up after returning from this function
//        string newPrefix = string(params.out_prefix) + "."  + NumberToString(MPIHelper::getInstance().getProcessID()) ;
//        params.out_prefix = (char *) newPrefix.c_str();
//    }

}

void usage(char* argv[]) {
    printCopyright(cout);
    cout << "Usage: " << argv[0] << " [OPTIONS] <file_name> [<output_file>]" << endl;
    cout << "GENERAL OPTIONS:" << endl;
    cout << "  -hh               Print this help dialog" << endl;
    cout << "  -h                Print help options for phylogenetic inference" << endl;
    cout << "  <file_name>       User tree in NEWICK format or split network in NEXUS format" << endl;
    cout << "  <output_file>     Output file to store results, default is '<file_name>.pda'" << endl;
    cout << "  -k <num_taxa>     Find optimal set of size <num_taxa>" << endl;
    cout << "  -k <min>:<max>    Find optimal sets of size from <min> to <max>" << endl;
    cout << "  -k <min>:<max>:<step>" << endl;
    cout << "                    Find optimal sets of size min, min+step, min+2*step,..." << endl;
    cout << "  -o <taxon>        Root name to compute rooted PD (default: unrooted)" << endl;
    cout << "  -if <file>        File containing taxa to be included into optimal sets" << endl;
    cout << "  -e <file>         File containing branch/split scale and taxa weights" << endl;
    cout << "  -all              Identify all multiple optimal sets" << endl;
    cout << "  -lim <max_limit>  The maximum number of optimal sets for each k if -a is specified" << endl;
    cout << "  -min              Compute minimal sets (default: maximal)" << endl;
    cout << "  -1out             Print taxa sets and scores to separate files" << endl;
    cout << "  -oldout           Print output compatible with version 0.3" << endl;
    cout << "  -v                Verbose mode" << endl;
    cout << endl;
    cout << "OPTIONS FOR PHYLOGENETIC DIVERSITY (PD):" << endl;
    cout << "  -root             Make the tree ROOTED, default is unrooted" << endl;
    cout << "    NOTE: this option and -o <taxon> cannot be both specified" << endl;
    cout << "  -g                Run greedy algorithm only (default: auto)" << endl;
    cout << "  -pr               Run pruning algorithm only (default: auto)" << endl;
    cout << endl;
    /*
    cout << "OPTIONS FOR SPLIT DIVERSITY:" << endl;
    cout << "  -exhaust          Force to use exhaustive search" << endl;
    cout << "    NOTE: by default, the program applies dynamic programming algorithm" << endl;
    cout << "          on circular networks and exhaustive search on general networks" << endl;
    cout << endl;*/
    cout << "OPTIONS FOR BUDGET CONSTRAINTS:" << endl;
    cout << "  -u <file>         File containing total budget and taxa preservation costs" << endl;
    cout << "  -b <budget>       Total budget to conserve taxa" << endl;
    cout << "  -b <min>:<max>    Find all sets with budget from <min> to <max>" << endl;
    cout << "  -b <min>:<max>:<step>" << endl;
    cout << "                    Find optimal sets with budget min, min+step, min+2*step,..." << endl;
    cout << endl;
    cout << "OPTIONS FOR AREA ANALYSIS:" << endl;
    cout << "  -ts <taxa_file>   Compute/maximize PD/SD of areas (combine with -k to maximize)" << endl;
    cout << "  -excl             Compute exclusive PD/SD" << endl;
    cout << "  -endem            Compute endemic PD/SD" << endl;
    cout << "  -compl <areas>    Compute complementary PD/SD given the listed <areas>" << endl;
    cout << endl;

    cout << "OPTIONS FOR VIABILITY CONSTRAINTS:" << endl;
    cout << "  -eco <food_web>   File containing food web matrix" << endl;
    cout << "  -k% <n>           Find optimal set of size relative the total number of taxa" << endl;
    cout << "  -diet <min_diet>  Minimum diet portion (%) to be preserved for each predator" << endl;
    cout << endl;
    //if (!full_command) exit(0);

    cout << "MISCELLANEOUS:" << endl;
    cout << "  -dd <sample_size> Compute PD distribution of random sets of size k" << endl;
    /*
    cout << "  -gbo <sitelh_file> Compute and output the alignment of (normalized)" << endl;
    cout << "                    expected frequencies given in site_ll_file" << endl;
	*/

    //	cout << "  -rep <times>        Repeat algorithm a number of times." << endl;
    //	cout << "  -noout              Print no output file." << endl;
    cout << endl;
    //cout << "HIDDEN OPTIONS: see the source code file pda.cpp::parseArg()" << endl;

    exit(0);
}

void usage_alisim(){
    cout << endl << "ALISIM: ALIGNMENT SIMULATOR" << endl
    << endl << "Usage: iqtree --alisim <OUTPUT_PREFIX> [-m MODEL] [-t TREE] ..." << endl << endl
    << "  --alisim OUTPUT_ALIGNMENT Activate AliSim and specify the output alignment filename"<< endl
    << "  -t TREE_FILE              Set the input tree file name" << endl
    << "  --length LENGTH           Set the length of the root sequence" << endl
    << "  --num-alignments NUMBER   Set the number of output datasets" << endl
    << "  --seqtype STRING          BIN, DNA, AA, CODON, MORPH{NUM_STATES} (default: auto-detect)" << endl
    << "                            For morphological data, 0<NUM_STATES<=32" << endl
    << "  --m MODEL_STRING          Specify the evolutionary model. See Manual for more detail" << endl
    << "  --mdef FILE               Name of a NEXUS model file to define new models (see Manual)" << endl
    << "  --fundi TAXA_LIST,RHO     Specify a list of taxa, and Rho (Fundi weight) for FunDi model" << endl
    << "  --indel <INS>,<DEL>       Set the insertion and deletion rate of the indel model,"<< endl
    << "                            relative to the substitution rate"<< endl
    << "  --indel-size <INS_DIS>,<DEL_DIS> Set the insertion and deletion size distributions" << endl
    << "  --sub-level-mixture       Enable the feature to simulate substitution-level mixture model"<< endl
    << "  --no-unaligned            Disable outputing a file of unaligned sequences "<< endl
    << "                            when using indel models"<< endl
    << "  --root-seq FILE,SEQ_NAME  Specify the root sequence from an alignment" << endl
    << "  -s FILE                   Specify the input sequence alignment" << endl
    << "  --no-copy-gaps            Disable copying gaps from input alignment (default: false)" << endl
    << "  --site-freq <OPTION>      Specify the option (MEAN (default), or SAMPLING, or MODEL)"<< endl
    << "                            to mimic the site-frequencies for mixture models from"<< endl
    << "                            the input alignment (see Manual)"<< endl
    << "  --site-rate <OPTION>      Specify the option (MEAN (default), or SAMPLING, or MODEL)"<< endl
    << "                            to mimic the discrete rate heterogeneity from"<< endl
    << "                            the input alignment (see Manual)"<< endl
    << "  -t RANDOM{MODEL,NUM_TAXA} Specify the model and the number of taxa to generate a random tree" << endl
    << "  -rlen MIN MEAN MAX        Specify three numbers: minimum, mean and maximum branch lengths" << endl
    << "                            when generating a random tree" << endl
    << "  -p FILE                   NEXUS/RAxML partition file" << endl
    << "                            Edge-linked proportional partition model" << endl
    << "  -q FILE                   Like -p but edge-linked equal partition model " << endl
    << "  -Q FILE                   Like -p but edge-unlinked partition model" << endl
    << "  --distribution FILE       Supply a definition file of distributions," << endl
    << "                            which could be used to generate random model parameters" << endl
    << "  --branch-distribution DIS Specify a distribution, from which branch lengths of the input trees" << endl
    << "                            are randomly generated and overridden." << endl
    << "  --branch-scale SCALE      Specify a value to scale all branch lengths" << endl
    << "  --single-output           Output all alignments into a single file" << endl
    << "  --write-all               Enable outputting internal sequences" << endl
    << "  --seed NUM                Random seed number (default: CPU clock)" << endl
    << "                            Be careful to make the AliSim reproducible," << endl
    << "                            users should specify the seed number" << endl
    << "  -gz                       Enable output compression but taking longer running time" << endl
    << "  -af phy|fasta             Set the output format (default: phylip)" << endl
    << "  User Manual is available at http://www.iqtree.org/doc/AliSim" << endl;
}

// TODO TD: add description for option use-nn-model
void usage_iqtree(char* argv[], bool full_command) {
    printCopyright(cout);
    cout << "Usage: iqtree [-s ALIGNMENT] [-p PARTITION] [-m MODEL] [-t TREE] ..." << endl << endl;
    cout << "GENERAL OPTIONS:" << endl
    << "  -h, --help           Print (more) help usages" << endl
    << "  -s FILE[,...,FILE]   PHYLIP/FASTA/NEXUS/CLUSTAL/MSF alignment file(s)" << endl
    << "  -s DIR               Directory of alignment files" << endl
    << "  --seqtype STRING     BIN, DNA, AA, NT2AA, CODON, MORPH (default: auto-detect)" << endl
    << "  -t FILE|PARS|RAND    Starting tree (default: 99 parsimony and BIONJ)" << endl
    << "  -o TAX[,...,TAX]     Outgroup taxon (list) for writing .treefile" << endl
    << "  --prefix STRING      Prefix for all output files (default: aln/partition)" << endl
    << "  --seed NUM           Random seed number, normally used for debugging purpose" << endl
    << "  --safe               Safe likelihood kernel to avoid numerical underflow" << endl
    << "  --mem NUM[G|M|%]     Maximal RAM usage in GB | MB | %" << endl
    << "  --runs NUM           Number of indepedent runs (default: 1)" << endl
    << "  -v, --verbose        Verbose mode, printing more messages to screen" << endl
    << "  -V, --version        Display version number" << endl
    << "  --quiet              Quiet mode, suppress printing to screen (stdout)" << endl
    << "  -fconst f1,...,fN    Add constant patterns into alignment (N=no. states)" << endl
    << "  --epsilon NUM        Likelihood epsilon for parameter estimate (default 0.01)" << endl
#ifdef _OPENMP
    << "  -T NUM|AUTO          No. cores/threads or AUTO-detect (default: 1)" << endl
    << "  --threads-max NUM    Max number of threads for -T AUTO (default: all cores)" << endl
#endif
    << endl << "CHECKPOINT:" << endl
    << "  --redo               Redo both ModelFinder and tree search" << endl
    << "  --redo-tree          Restore ModelFinder and only redo tree search" << endl
    << "  --undo               Revoke finished run, used when changing some options" << endl
    << "  --cptime NUM         Minimum checkpoint interval (default: 60 sec and adapt)" << endl
    << endl << "PARTITION MODEL:" << endl
    << "  -p FILE|DIR          NEXUS/RAxML partition file or directory with alignments" << endl
    << "                       Edge-linked proportional partition model" << endl
    << "  -q FILE|DIR          Like -p but edge-linked equal partition model " << endl
    << "  -Q FILE|DIR          Like -p but edge-unlinked partition model" << endl
    << "  -S FILE|DIR          Like -p but separate tree inference" << endl
    << "  --subsample NUM      Randomly sub-sample partitions (negative for complement)" << endl
    << "  --subsample-seed NUM Random number seed for --subsample" << endl
    << endl << "LIKELIHOOD/QUARTET MAPPING:" << endl
    << "  --lmap NUM           Number of quartets for likelihood mapping analysis" << endl
    << "  --lmclust FILE       NEXUS file containing clusters for likelihood mapping" << endl
    << "  --quartetlh          Print quartet log-likelihoods to .quartetlh file" << endl
    << endl << "TREE SEARCH ALGORITHM:" << endl
//            << "  -pll                 Use phylogenetic likelihood library (PLL) (default: off)" << endl
    << "  --ninit NUM          Number of initial parsimony trees (default: 100)" << endl
    << "  --ntop NUM           Number of top initial trees (default: 20)" << endl
    << "  --nbest NUM          Number of best trees retained during search (defaut: 5)" << endl
    << "  -n NUM               Fix number of iterations to stop (default: OFF)" << endl
    << "  --nstop NUM          Number of unsuccessful iterations to stop (default: 100)" << endl
    << "  --perturb NUM        Perturbation strength for randomized NNI (default: 0.5)" << endl
    << "  --radius NUM         Radius for parsimony SPR search (default: 6)" << endl
    << "  --allnni             Perform more thorough NNI search (default: OFF)" << endl
    << "  -g FILE              (Multifurcating) topological constraint tree file" << endl
    << "  --fast               Fast search to resemble FastTree" << endl
    << "  --polytomy           Collapse near-zero branches into polytomy" << endl
    << "  --tree-fix           Fix -t tree (no tree search performed)" << endl
    << "  --treels             Write locally optimal trees into .treels file" << endl
    << "  --show-lh            Compute tree likelihood without optimisation" << endl
#ifdef IQTREE_TERRAPHAST
    << "  --terrace            Check if the tree lies on a phylogenetic terrace" << endl
#endif
//            << "  -iqp                 Use the IQP tree perturbation (default: randomized NNI)" << endl
//            << "  -iqpnni              Switch back to the old IQPNNI tree search algorithm" << endl
    << endl << "ULTRAFAST BOOTSTRAP/JACKKNIFE:" << endl
    << "  -B, --ufboot NUM     Replicates for ultrafast bootstrap (>=1000)" << endl
    << "  -J, --ufjack NUM     Replicates for ultrafast jackknife (>=1000)" << endl
    << "  --jack-prop NUM      Subsampling proportion for jackknife (default: 0.5)" << endl
    << "  --sampling STRING    GENE|GENESITE resampling for partitions (default: SITE)" << endl
    << "  --boot-trees         Write bootstrap trees to .ufboot file (default: none)" << endl
    << "  --wbtl               Like --boot-trees but also writing branch lengths" << endl
//            << "  -n <#iterations>     Minimum number of iterations (default: 100)" << endl
    << "  --nmax NUM           Maximum number of iterations (default: 1000)" << endl
    << "  --nstep NUM          Iterations for UFBoot stopping rule (default: 100)" << endl
    << "  --bcor NUM           Minimum correlation coefficient (default: 0.99)" << endl
    << "  --beps NUM           RELL epsilon to break tie (default: 0.5)" << endl
    << "  --bnni               Optimize UFBoot trees by NNI on bootstrap alignment" << endl
    << endl << "NON-PARAMETRIC BOOTSTRAP/JACKKNIFE:" << endl
    << "  -b, --boot NUM       Replicates for bootstrap + ML tree + consensus tree" << endl
    << "  -j, --jack NUM       Replicates for jackknife + ML tree + consensus tree" << endl
    << "  --jack-prop NUM      Subsampling proportion for jackknife (default: 0.5)" << endl
    << "  --bcon NUM           Replicates for bootstrap + consensus tree" << endl
    << "  --bonly NUM          Replicates for bootstrap only" << endl
#ifdef USE_BOOSTER
    << "  --tbe                Transfer bootstrap expectation" << endl
#endif
//            << "  -t <threshold>       Minimum bootstrap support [0...1) for consensus tree" << endl
    << endl << "SINGLE BRANCH TEST:" << endl
    << "  --alrt NUM           Replicates for SH approximate likelihood ratio test" << endl
    << "  --alrt 0             Parametric aLRT test (Anisimova and Gascuel 2006)" << endl
    << "  --abayes             approximate Bayes test (Anisimova et al. 2011)" << endl
    << "  --lbp NUM            Replicates for fast local bootstrap probabilities" << endl
    << endl << "MODEL-FINDER:" << endl
    << "  --use-nn-model       Use neural network for tree inference" << endl
    << "  --nn-path-model      Neural network file for substitution model (onnx format)" << endl
    << "  --nn-path-rates      Neural network file for alpha value (onnx format)" << endl
    << "  -m TESTONLY          Standard model selection (like jModelTest, ProtTest)" << endl
    << "  -m TEST              Standard model selection followed by tree inference" << endl
    << "  -m MF                Extended model selection with FreeRate heterogeneity" << endl
    << "  -m MFP               Extended model selection followed by tree inference" << endl
    << "  -m ...+LM            Additionally test Lie Markov models" << endl
    << "  -m ...+LMRY          Additionally test Lie Markov models with RY symmetry" << endl
    << "  -m ...+LMWS          Additionally test Lie Markov models with WS symmetry" << endl
    << "  -m ...+LMMK          Additionally test Lie Markov models with MK symmetry" << endl
    << "  -m ...+LMSS          Additionally test strand-symmetric models" << endl
    << "  --mset STRING        Restrict search to models supported by other programs" << endl
    << "                       (raxml, phyml, mrbayes, beast1 or beast2)" << endl
    << "  --mset STR,...       Comma-separated model list (e.g. -mset WAG,LG,JTT)" << endl
    << "  --msub STRING        Amino-acid model source" << endl
    << "                       (nuclear, mitochondrial, chloroplast or viral)" << endl
    << "  --mfreq STR,...      List of state frequencies" << endl
    << "  --mrate STR,...      List of rate heterogeneity among sites" << endl
    << "                       (e.g. -mrate E,I,G,I+G,R is used for -m MF)" << endl
    << "  --cmin NUM           Min categories for FreeRate model [+R] (default: 2)" << endl
    << "  --cmax NUM           Max categories for FreeRate model [+R] (default: 10)" << endl
    << "  --merit AIC|AICc|BIC  Akaike|Bayesian information criterion (default: BIC)" << endl
//            << "  -msep                Perform model selection and then rate selection" << endl
    << "  --mtree              Perform full tree search for every model" << endl
    << "  --madd STR,...       List of mixture models to consider" << endl
    << "  --mdef FILE          Model definition NEXUS file (see Manual)" << endl
    << "  --modelomatic        Find best codon/protein/DNA models (Whelan et al. 2015)" << endl

    << endl << "PARTITION-FINDER:" << endl
    << "  --merge              Merge partitions to increase model fit" << endl
    << "  --merge greedy|rcluster|rclusterf" << endl
    << "                       Set merging algorithm (default: rclusterf)" << endl
    << "  --merge-model 1|all  Use only 1 or all models for merging (default: 1)" << endl
    << "  --merge-model STR,..." << endl
    << "                       Comma-separated model list for merging" << endl
    << "  --merge-rate 1|all   Use only 1 or all rate heterogeneity (default: 1)" << endl
    << "  --merge-rate STR,..." << endl
    << "                       Comma-separated rate list for merging" << endl
    << "  --rcluster NUM       Percentage of partition pairs for rcluster algorithm" << endl
    << "  --rclusterf NUM      Percentage of partition pairs for rclusterf algorithm" << endl
    << "  --rcluster-max NUM   Max number of partition pairs (default: 10*partitions)" << endl

    << endl << "SUBSTITUTION MODEL:" << endl
    << "  -m STRING            Model name string (e.g. GTR+F+I+G)" << endl
    << "                 DNA:  HKY (default), JC, F81, K2P, K3P, K81uf, TN/TrN, TNef," << endl
    << "                       TIM, TIMef, TVM, TVMef, SYM, GTR, or 6-digit model" << endl
    << "                       specification (e.g., 010010 = HKY)" << endl
    << "             Protein:  LG (default), Poisson, cpREV, mtREV, Dayhoff, mtMAM," << endl
    << "                       JTT, WAG, mtART, mtZOA, VT, rtREV, DCMut, PMB, HIVb," << endl
    << "                       HIVw, JTTDCMut, FLU, Blosum62, GTR20, mtMet, mtVer, mtInv, FLAVI," << endl
    << "			Q.LG, Q.pfam, Q.pfam_gb, Q.bird, Q.mammal, Q.insect, Q.plant, Q.yeast" << endl
    << "     Protein mixture:  C10,...,C60, EX2, EX3, EHO, UL2, UL3, EX_EHO, LG4M, LG4X" << endl
    << "              Binary:  JC2 (default), GTR2" << endl
    << "     Empirical codon:  KOSI07, SCHN05" << endl
    << "   Mechanistic codon:  GY (default), MG, MGK, GY0K, GY1KTS, GY1KTV, GY2K," << endl
    << "                       MG1KTS, MG1KTV, MG2K" << endl
    << "Semi-empirical codon:  XX_YY where XX is empirical and YY is mechanistic model" << endl
    << "      Morphology/SNP:  MK (default), ORDERED, GTR" << endl
    << "      Lie Markov DNA:  1.1, 2.2b, 3.3a, 3.3b, 3.3c, 3.4, 4.4a, 4.4b, 4.5a," << endl
    << "                       4.5b, 5.6a, 5.6b, 5.7a, 5.7b, 5.7c, 5.11a, 5.11b, 5.11c," << endl
    << "                       5.16, 6.6, 6.7a, 6.7b, 6.8a, 6.8b, 6.17a, 6.17b, 8.8," << endl
    << "                       8.10a, 8.10b, 8.16, 8.17, 8.18, 9.20a, 9.20b, 10.12," << endl
    << "                       10.34, 12.12 (optionally prefixed by RY, WS or MK)" << endl
    << "      Non-reversible:  STRSYM (strand symmetric model, equiv. WS6.6)," << endl
    << "                       NONREV, UNREST (unrestricted model, equiv. 12.12)" << endl
    << "                       NQ.pfam, NQ.bird, NQ.mammal, NQ.insect, NQ.plant, NQ.yeast" << endl
    << "           Otherwise:  Name of file containing user-model parameters" << endl
    << endl << "STATE FREQUENCY:" << endl
    << "  -m ...+F             Empirically counted frequencies from alignment" << endl
    << "  -m ...+FO            Optimized frequencies by maximum-likelihood" << endl
    << "  -m ...+FQ            Equal frequencies" << endl
    << "  -m ...+FRY           For DNA, freq(A+G)=1/2=freq(C+T)" << endl
    << "  -m ...+FWS           For DNA, freq(A+T)=1/2=freq(C+G)" << endl
    << "  -m ...+FMK           For DNA, freq(A+C)=1/2=freq(G+T)" << endl
    << "  -m ...+Fabcd         4-digit constraint on ACGT frequency" << endl
    << "                       (e.g. +F1221 means f_A=f_T, f_C=f_G)" << endl
    << "  -m ...+FU            Amino-acid frequencies given protein matrix" << endl
    << "  -m ...+F1x4          Equal NT frequencies over three codon positions" << endl
    << "  -m ...+F3x4          Unequal NT frequencies over three codon positions" << endl

    << endl << "RATE HETEROGENEITY AMONG SITES:" << endl
    << "  -m ...+I             A proportion of invariable sites" << endl
    << "  -m ...+G[n]          Discrete Gamma model with n categories (default n=4)" << endl
    << "  -m ...*G[n]          Discrete Gamma model with unlinked model parameters" << endl
    << "  -m ...+I+G[n]        Invariable sites plus Gamma model with n categories" << endl
    << "  -m ...+R[n]          FreeRate model with n categories (default n=4)" << endl
    << "  -m ...*R[n]          FreeRate model with unlinked model parameters" << endl
    << "  -m ...+I+R[n]        Invariable sites plus FreeRate model with n categories" << endl
    << "  -m ...+Hn            Heterotachy model with n classes" << endl
    << "  -m ...*Hn            Heterotachy model with n classes and unlinked parameters" << endl
    << "  --alpha-min NUM      Min Gamma shape parameter for site rates (default: 0.02)" << endl
    << "  --gamma-median       Median approximation for +G site rates (default: mean)" << endl
    << "  --rate               Write empirical Bayesian site rates to .rate file" << endl
    << "  --mlrate             Write maximum likelihood site rates to .mlrate file" << endl
//            << "  --mhrate             Computing site-specific rates to .mhrate file using" << endl
//            << "                       Meyer & von Haeseler (2003) method" << endl

    << endl << "POLYMORPHISM AWARE MODELS (PoMo):"                                           << endl
    << "  -s FILE              Input counts file (see manual)"                               << endl
    << "  -m ...+P             DNA substitution model (see above) used with PoMo"            << endl
    << "  -m ...+N<POPSIZE>    Virtual population size (default: 9)"                         << endl
// TODO DS: Maybe change default to +WH.
    << "  -m ...+WB|WH|S]      Weighted binomial sampling"       << endl
    << "  -m ...+WH            Weighted hypergeometric sampling" << endl
    << "  -m ...+S             Sampled sampling"              << endl
    << "  -m ...+G[n]          Discrete Gamma rate with n categories (default n=4)"    << endl
// TODO DS: Maybe change default to +WH.

    << endl << "COMPLEX MODELS:" << endl
    << "  -m \"MIX{m1,...,mK}\"  Mixture model with K components" << endl
    << "  -m \"FMIX{f1,...fK}\"  Frequency mixture model with K components" << endl
    << "  --mix-opt            Optimize mixture weights (default: detect)" << endl
    << "  -m ...+ASC           Ascertainment bias correction" << endl
    << "  --tree-freq FILE     Input tree to infer site frequency model" << endl
    << "  --site-freq FILE     Input site frequency model file" << endl
    << "  --freq-max           Posterior maximum instead of mean approximation" << endl

    << endl << "TREE TOPOLOGY TEST:" << endl
    << "  --trees FILE         Set of trees to evaluate log-likelihoods" << endl
    << "  --test NUM           Replicates for topology test" << endl
    << "  --test-weight        Perform weighted KH and SH tests" << endl
    << "  --test-au            Approximately unbiased (AU) test (Shimodaira 2002)" << endl
    << "  --sitelh             Write site log-likelihoods to .sitelh file" << endl

    << endl << "ANCESTRAL STATE RECONSTRUCTION:" << endl
    << "  --ancestral          Ancestral state reconstruction by empirical Bayes" << endl
    << "  --asr-min NUM        Min probability of ancestral state (default: equil freq)" << endl

    << endl << "TEST OF SYMMETRY:" << endl
    << "  --symtest               Perform three tests of symmetry" << endl
    << "  --symtest-only          Do --symtest then exist" << endl
//    << "  --bisymtest             Perform three binomial tests of symmetry" << endl
//    << "  --symtest-perm NUM      Replicates for permutation tests of symmetry" << endl
    << "  --symtest-remove-bad    Do --symtest and remove bad partitions" << endl
    << "  --symtest-remove-good   Do --symtest and remove good partitions" << endl
    << "  --symtest-type MAR|INT  Use MARginal/INTernal test when removing partitions" << endl
    << "  --symtest-pval NUMER    P-value cutoff (default: 0.05)" << endl
    << "  --symtest-keep-zero     Keep NAs in the tests" << endl

    << endl << "CONCORDANCE FACTOR ANALYSIS:" << endl
    << "  -t FILE              Reference tree to assign concordance factor" << endl
    << "  --gcf FILE           Set of source trees for gene concordance factor (gCF)" << endl
    << "  --df-tree            Write discordant trees associated with gDF1" << endl
    << "  --scf NUM            Number of quartets for site concordance factor (sCF)" << endl
    << "  --scfl NUM           Like --scf but using likelihood (recommended)" << endl
    << "  -s FILE              Sequence alignment for --scf" << endl
    << "  -p FILE|DIR          Partition file or directory for --scf" << endl
    << "  --cf-verbose         Write CF per tree/locus to cf.stat_tree/_loci" << endl
    << "  --cf-quartet         Write sCF for all resampled quartets to .cf.quartet" << endl;
    
    usage_alisim();
    cout

    << endl << "ANALYSIS WITH GENTRIUS ALGORITHM:" << endl
    << "  --gentrius FILE      File must contain either a single species-tree or a set of subtrees." << endl
    << "  -pr_ab_matrix FILE   Presence-absence matrix of loci coverage." << endl
    << "  -s FILE              PHYLIP/FASTA/NEXUS/CLUSTAL/MSF alignment file(s)" << endl
    << "  -p FILE              NEXUS/RAxML partition file" << endl
    << "  -g_stop_t NUM        Stop after NUM species-trees were generated, or use 0 to turn off this stopping rule. Default: 1MLN trees."<< endl
    << "  -g_stop_i NUM        Stop after NUM intermediate trees were visited, or use 0 to turn off this stopping rule. Default: 10MLN trees." << endl
    << "  -g_stop_h NUM        Stop after NUM hours (CPU time), or use 0 to turn off this stopping rule. Default: 7 days." << endl
    << "  -g_non_stop          Turn off all stopping rules." << endl
    << "  -g_query FILE        Species-trees to test for identical set of subtrees." << endl
    << "  -g_print             Write all generated species-trees. WARNING: there might be millions of trees!" << endl
    << "  -g_print_lim NUM     Limit on the number of species-trees to be written." << endl
    << "  -g_print_induced     Write induced partition subtrees." << endl
    << "  -g_print_m           Write presence-absence matrix." << endl
    << "  -g_rm_leaves NUM     Invoke reverse analysis for complex datasets." << endl
    
    << endl << "GENOMIC EPIDEMIOLOGICAL ANALYSIS:" << endl
    << "  --pathogen           Apply CMAPLE tree search algorithm if sequence" << endl
    << "                       divergence is low, otherwise, apply IQ-TREE algorithm." << endl
    << "  --pathogen-force     Apply CMAPLE tree search algorithm regardless" << endl
    << "                       of sequence divergence." << endl
    << "  --alrt <num_rep>     Specify number of replicates to compute SH-aLRT." << endl
    << "  --sprta              Compute SPRTA (DeMaio et al., 2024) branch supports." << endl
    << "  --sprta-zero-branch  Compute SPRTA supports for zero-length branches." << endl
    << "  --sprta-other-places Output alternative SPRs and their SPRTA supports." << endl
    << "  -T <num_thread>      Specify number of threads used for computing" << endl
    << "                       branch supports (SH-aLRT or SPRTA)." << endl


    
#ifdef USE_LSD2
    << endl << "TIME TREE RECONSTRUCTION:" << endl
    << "  --date FILE          File containing dates of tips or ancestral nodes" << endl
    << "  --date TAXNAME       Extract dates from taxon names after last '|'" << endl
    << "  --date-tip STRING    Tip dates as a real number or YYYY-MM-DD" << endl
    << "  --date-root STRING   Root date as a real number or YYYY-MM-DD" << endl
    << "  --date-ci NUM        Number of replicates to compute confidence interval" << endl
    << "  --clock-sd NUM       Std-dev for lognormal relaxed clock (default: 0.2)" << endl
    << "  --date-no-outgroup   Exclude outgroup from time tree" << endl
    << "  --date-outlier NUM   Z-score cutoff to remove outlier tips/nodes (e.g. 3)" << endl
    << "  --date-options \"..\"  Extra options passing directly to LSD2" << endl
    << "  --dating STRING      Dating method: LSD for least square dating (default)" << endl
#endif
    << endl;
    

//            << endl << "TEST OF MODEL HOMOGENEITY:" << endl
//            << "  -m WHTEST            Testing model (GTR+G) homogeneity assumption using" << endl
//            << "                       Weiss & von Haeseler (2003) method" << endl
//            << "  -ns <#simulations>   #Simulations to obtain null-distribution (default: 1000)" << endl

    if (full_command)
    cout
        << endl << "CONSENSUS RECONSTRUCTION:" << endl
        << "  -t FILE              Set of input trees for consensus reconstruction" << endl
        << "  --sup-min NUM        Min split support, 0.5 for majority-rule consensus" << endl
        << "                       (default: 0, extended consensus)" << endl
        << "  --burnin NUM         Burnin number of trees to ignore" << endl
        << "  --con-tree           Compute consensus tree to .contree file" << endl
        << "  --con-net            Computing consensus network to .nex file" << endl
        << "  --support FILE       Assign support values into this tree from -t trees" << endl
        //<< "  -sup2 FILE           Like -sup but -t trees can have unequal taxon sets" << endl
        << "  --suptag STRING      Node name (or ALL) to assign tree IDs where node occurs" << endl
        << endl << "TREE DISTANCE BY ROBINSON-FOULDS (RF) METRIC:" << endl
        << "  --tree-dist-all      Compute all-to-all RF distances for -t trees" << endl
        << "  --tree-dist FILE     Compute RF distances between -t trees and this set" << endl
        << "  --tree-dist2 FILE    Like -rf but trees can have unequal taxon sets" << endl
    //            << "  -rf_adj              Computing RF distances of adjacent trees in <treefile>" << endl
    //            << "  -wja                 Write ancestral sequences by joint reconstruction" << endl


        << endl

        << "GENERATING RANDOM TREES:" << endl
        << "  -r NUM               No. taxa for Yule-Harding random tree" << endl
        << "  --rand UNI|CAT|BAL   UNIform | CATerpillar | BALanced random tree" << endl
        //<< "  --rand NET           Random circular split network" << endl
        << "  --rlen NUM NUM NUM   min, mean, and max random branch lengths" << endl

        << endl << "MISCELLANEOUS:" << endl
        << "  --keep-ident         Keep identical sequences (default: remove & finally add)" << endl
        << "  -blfix               Fix branch lengths of user tree passed via -te" << endl
        << "  -blscale             Scale branch lengths of user tree passed via -t" << endl
        << "  -blmin               Min branch length for optimization (default 0.000001)" << endl
        << "  -blmax               Max branch length for optimization (default 100)" << endl
        << "  -wslr                Write site log-likelihoods per rate category" << endl
        << "  -wslm                Write site log-likelihoods per mixture class" << endl
        << "  -wslmr               Write site log-likelihoods per mixture+rate class" << endl
        << "  -wspr                Write site probabilities per rate category" << endl
        << "  -wspm                Write site probabilities per mixture class" << endl
        << "  -wspmr               Write site probabilities per mixture+rate class" << endl
        << "  --partlh             Write partition log-likelihoods to .partlh file" << endl
        << "  --no-outfiles        Suppress printing output files" << endl
        << "  --eigenlib           Use Eigen3 library" << endl
        << "  -alninfo             Print alignment sites statistics to .alninfo" << endl
    //            << "  -d <file>            Reading genetic distances from file (default: JC)" << endl
    //			<< "  -d <outfile>         Calculate the distance matrix inferred from tree" << endl
    //			<< "  -stats <outfile>     Output some statistics about branch lengths" << endl
    //			<< "  -comp <treefile>     Compare tree with each in the input trees" << endl;


        << endl;

    if (full_command) {
        //TODO Print other options here (to be added)
    }

    exit(0);
}

void quickStartGuide() {
    printCopyright(cout);
    cout << "Command-line examples (replace 'iqtree3 ...' by actual path to executable):" << endl << endl
         << "1. Infer maximum-likelihood tree from a sequence alignment (example.phy)" << endl
         << "   with the best-fit model automatically selected by ModelFinder:" << endl
         << "     iqtree3 -s example.phy" << endl << endl
         << "2. Perform ModelFinder without subsequent tree inference:" << endl
         << "     iqtree3 -s example.phy -m MF" << endl
         << "   (use '-m TEST' to resemble jModelTest/ProtTest)" << endl << endl
         << "3. Combine ModelFinder, tree search, ultrafast bootstrap and SH-aLRT test:" << endl
         << "     iqtree3 -s example.phy --alrt 1000 -B 1000" << endl << endl
         << "4. Perform edge-linked proportional partition model (example.nex):" << endl
         << "     iqtree3 -s example.phy -p example.nex" << endl
         << "   (replace '-p' by '-Q' for edge-unlinked model)" << endl << endl
         << "5. Find best partition scheme by possibly merging partitions:" << endl
         << "     iqtree3 -s example.phy -p example.nex -m MF+MERGE" << endl
         << "   (use '-m TESTMERGEONLY' to resemble PartitionFinder)" << endl << endl
         << "6. Find best partition scheme followed by tree inference and bootstrap:" << endl
         << "     iqtree3 -s example.phy -p example.nex -m MFP+MERGE -B 1000" << endl << endl
#ifdef _OPENMP
         << "7. Use 4 CPU cores to speed up computation: add '-T 4' option" << endl << endl
#endif
         << "8. Polymorphism-aware model with HKY nucleotide model and Gamma rate:" << endl
         << "     iqtree3 -s counts_file.cf -m HKY+P+G" << endl << endl
    // BQM: this example is too complicated
//         << "9. PoMo mixture with virtual popsize 5 and weighted binomial sampling:" << endl
//         << "     iqtree3 -s counts_file.cf -m \"MIX{HKY+P{EMP},JC+P}+N5+WB\"" << endl << endl
         << "To show all available options: run 'iqtree3 -h'" << endl << endl
         << "Have a look at the tutorial and manual for more information:" << endl
         << "     http://www.iqtree.org" << endl << endl;
    exit(0);
}

InputType detectInputFile(const char *input_file) {

    if (!fileExists(input_file))
        outError("File not found ", input_file);

    try {
        igzstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(input_file);

        unsigned char ch, ch2;
        int count = 0;
        do {
            in >> ch;
        } while (ch <= 32 && !in.eof() && count++ < 20);
        in >> ch2;
        in.close();
        switch (ch) {
            case '#': return IN_NEXUS;
            case '(': return IN_NEWICK;
            case '[': return IN_NEWICK;
            case '>': return IN_FASTA;
            case 'C': if (ch2 == 'L') return IN_CLUSTAL;
                      else if (ch2 == 'O') return IN_COUNTS;
                      else return IN_OTHER;
            case '!': if (ch2 == '!') return IN_MSF; else return IN_OTHER;
            default:
                if (isdigit(ch)) return IN_PHYLIP;
                return IN_OTHER;
        }
    } catch (ios::failure) {
        outError("Cannot read file ", input_file);
    } catch (...) {
        outError("Cannot read file ", input_file);
    }
    return IN_OTHER;
}

bool overwriteFile(char *filename) {
    ifstream infile(filename);
    if (infile.is_open()) {
        cout << "Overwrite " << filename << " (y/n)? ";
        char ch;
        cin >> ch;
        if (ch != 'Y' && ch != 'y') {
            infile.close();
            return false;
        }
    }
    infile.close();
    return true;
}

void parseAreaName(char *area_names, set<string> &areas) {
    string all = area_names;
    int pos;
    while (!all.empty()) {
        pos = all.find(',');
        if (pos < 0) pos = all.length();
        areas.insert(all.substr(0, pos));
        if (pos >= (signed int) all.length())
            all = "";
        else
            all = all.substr(pos + 1);
    }
}

double logFac(const int num) {
    if (num < 0) return -1.0;
    if (num == 0) return 0.0;
    double ret = 0;
    for (int i = 1; i <= num; i++)
        ret += log((double) i);
    return ret;
}

template <typename I>
I random_element(I begin, I end)
{
    const unsigned long n = std::distance(begin, end);
    const unsigned long divisor = (RAND_MAX + 1) / n;

    unsigned long k;
    do { k = std::rand() / divisor; } while (k >= n);

    return std::advance(begin, k);
}

template <class T>
inline T quantile(const vector<T>& v, const double q) {
    unsigned int size = v.size();
    if (q <= 0) return *std::min_element(v.begin(), v.end());
    if (q >= 1) return *std::max_element(v.begin(), v.end());
    //double pos = (size - 1) * q;
    //unsigned int ind = (unsigned int)(pos);
    //double delta = pos - ind;
    vector<T> w(size);
    std::copy(v, v.begin() + size, w.begin());
}

#define RAN_STANDARD 1
#define RAN_SPRNG    2
#define RAN_RAND4    3

#define RAN_TYPE 2

#if RAN_TYPE == RAN_STANDARD

int init_random(int seed) {
    srand(seed);
    cout << "(Using rand() - Standard Random Number Generator)" << endl;
    return seed;
}

int finish_random() {
	return 0;
}


#elif RAN_TYPE == RAN_RAND4
/******************************************************************************/
/* random numbers generator  (Numerical recipes)                              */
/******************************************************************************/

/* variable */
long _idum;

/* definitions */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double randomunitintervall()
/* Long period (> 2e18) random number generator. Returns a uniform random
   deviate between 0.0 and 1.0 (exclusive of endpoint values).

   Source:
   Press et al., "Numerical recipes in C", Cambridge University Press, 1992
   (chapter 7 "Random numbers", ran2 random number generator) */ {
    int j;
    long k;
    static long _idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (_idum <= 0) {
        if (-(_idum) < 1)
            _idum = 1;
        else
            _idum = -(_idum);
        _idum2 = (_idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (_idum) / IQ1;
            _idum = IA1 * (_idum - k * IQ1) - k*IR1;
            if (_idum < 0)
                _idum += IM1;
            if (j < NTAB)
                iv[j] = _idum;
        }
        iy = iv[0];
    }
    k = (_idum) / IQ1;
    _idum = IA1 * (_idum - k * IQ1) - k*IR1;
    if (_idum < 0)
        _idum += IM1;
    k = _idum2 / IQ2;
    _idum2 = IA2 * (_idum2 - k * IQ2) - k*IR2;
    if (_idum2 < 0)
        _idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - _idum2;
    iv[j] = _idum;
    if (iy < 1)
        iy += IMM1;
    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
} /* randomunitintervall */

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

int init_random(int seed) /* RAND4 */ {
    //    srand((unsigned) time(nullptr));
    //    if (seed < 0)
    // 	seed = rand();
    _idum = -(long) seed;
#ifndef PARALLEL
    cout << "(Using RAND4 Random Number Generator)" << endl;
#else /* PARALLEL */
    {
        int n;
        if (PP_IamMaster) {
            cout << "(Using RAND4 Random Number Generator with leapfrog method)" << endl;
        }
        for (n = 0; n < PP_Myid; n++)
            (void) randomunitintervall();
        if (verbose_mode >= VB_MED) {
            cout << "(" << PP_Myid << ") !!! random seed set to " << seed << ", " << n << " drawn !!!" << endl;
        }
    }
#endif
    return (seed);
} /* initrandom */

int finish_random() {
	return 0;
}
/******************/

#else /* SPRNG */

/******************/

int *randstream;
/**
   vector of random streams for multiple threads
 **/
vector<int*> rstream_vec;
/**
   vector of generators for multimodal discrete distributions and continuous Gamma distributions
   Note that multimodal discrete distributions were used to randomly select a site
   when handling Indel/Sub events with the Gillespie algorithm.
 **/
vector<default_random_engine> generator_vec;

int init_random(int seed, bool write_info, int** rstream) {
    //    srand((unsigned) time(nullptr));
    if (seed < 0)
        seed = make_sprng_seed();
#ifndef PARALLEL
    if (write_info)
    	cout << "(Using SPRNG - Scalable Parallel Random Number Generator)" << endl;
    if (rstream) {
        *rstream = init_sprng(0, 1, seed, SPRNG_DEFAULT); /*init stream*/
    } else {
        randstream = init_sprng(0, 1, seed, SPRNG_DEFAULT); /*init stream*/
        if (verbose_mode >= VB_MED) {
            print_sprng(randstream);
        }
    }
#else /* PARALLEL */
    if (PP_IamMaster && write_info) {
        cout << "(Using SPRNG - Scalable Parallel Random Number Generator)" << endl;
    }
    /* MPI_Bcast(&seed, 1, MPI_UNSIGNED, PP_MyMaster, MPI_COMM_WORLD); */
    if (rstream) {
        *rstream = init_sprng(PP_Myid, PP_NumProcs, seed, SPRNG_DEFAULT); /*initialize stream*/
    } else {
        randstream = init_sprng(PP_Myid, PP_NumProcs, seed, SPRNG_DEFAULT); /*initialize stream*/
        if (verbose_mode >= VB_MED) {
            cout << "(" << PP_Myid << ") !!! random seed set to " << seed << " !!!" << endl;
            print_sprng(randstream);
        }
    }
#endif /* PARALLEL */
    return (seed);
} /* initrandom */

int finish_random(int *rstream) {
    if (rstream)
        return free_sprng(rstream);
    else
        return free_sprng(randstream);
}

int init_multi_rstreams()
{
#if RAN_TYPE == RAN_SPRNG
    #ifdef _OPENMP
    // get the number of all threads (not just physical)
    const int num_threads = countPhysicalCPUCores();
    
    // initialize a vector of random seeds
    size_t ran_seed_vec_size = 2*num_threads;
    size_t i;
    vector<int64_t> ran_seed_vec(ran_seed_vec_size);
    
    // generate a vector of random seeds
//    const int MAX_INT32 = 2147483647;
//    for (i = 0; i < ran_seed_vec_size; ++i)
//        ran_seed_vec[i] = random_int(MAX_INT32);
      for (i = 0; i < ran_seed_vec_size; ++i)
          ran_seed_vec[i] = Params::getInstance().ran_seed + (int64_t)MPIHelper::getInstance().getProcessID()*1000 + i+1;

    // print the random seeds for debugging purposes if needed
    if (verbose_mode >= VB_DEBUG)
    {
        std::cout << "- Random seeds used:";
        for (i = 0; i < ran_seed_vec_size; ++i)
            std::cout << " " << ran_seed_vec[i];
        std::cout << std::endl;
    }
        
    // initialize the vector of random streams
    rstream_vec.resize(num_threads);
    for (i = 0; i < num_threads; ++i) {
        init_random(ran_seed_vec[i], false, &rstream_vec[i]);
    }
        
    // initialize the continuous Gamma generators
    generator_vec.resize(num_threads);
    for (i = 0; i < num_threads; ++i)
        generator_vec[i].seed(ran_seed_vec[i+num_threads]);
    
    #endif
#else
    outError("Oops! Only SPRNG is now supported for the parallel version.")
#endif
    
    return rstream_vec.size();
}

int finish_multi_rstreams()
{
    // finish the random streams one by one
    for (size_t i = 0; i < rstream_vec.size(); ++i)
        finish_random(rstream_vec[i]);
    
    // clear the vector of random streams
    rstream_vec.clear();
    
    return rstream_vec.size();
}

#endif /* USE_SPRNG */

/******************/

/* returns a random integer in the range [0; n - 1] */
int random_int(int n, int *rstream) {
    return (int) floor(random_double(rstream) * n);
} /* randominteger */

/* returns a random integer in the range [a; b] */
int random_int(int a, int b) {
	ASSERT(b > a);
	//return a + (RAND_MAX * rand() + rand()) % (b + 1 - a);
	return a + random_int(b - a);
}

double random_double(int *rstream) {
#ifndef FIXEDINTRAND
#ifndef PARALLEL
#if RAN_TYPE == RAN_STANDARD
    return ((double) rand()) / ((double) RAND_MAX + 1);
#elif RAN_TYPE == RAN_SPRNG
    if (rstream)
        return sprng(rstream);
    else
        return sprng(randstream);
#else /* NO_SPRNG */
    return randomunitintervall();
#endif /* NO_SPRNG */
#else /* NOT PARALLEL */
#if RAN_TYPE == RAN_SPRNG
    if (rstream)
        return sprng(rstream);
    else
        return sprng(randstream);
#else /* NO_SPRNG */
    int m;
    for (m = 1; m < PP_NumProcs; m++)
        (void) randomunitintervall();
    PP_randn += (m - 1);
    PP_rand++;
    return randomunitintervall();
#endif /* NO_SPRNG */
#endif /* NOT PARALLEL */
#else /* FIXEDINTRAND */
    cerr << "!!! fixed \"random\" integers for testing purposes !!!" << endl;
    return 0.0;
#endif /* FIXEDINTRAND */

}

/**
 * returns a random double based on an exponential distribution
 * @param mean the mean of exponential distribution
 */
double random_double_exponential_distribution(double mean, int *rstream)
{
    double ran;
    do
        ran = random_double(rstream);
    while (ran == 0.0);
    return -mean * log(ran);
}

/**
 * geometric random number generation (starting from 0)
 * Modified from W. Fletcher and Z. Yang, “INDELible: A flexible simulator of biological sequence evolution,” Mol. Biol. Evol., vol. 26, no. 8, pp. 1879–1888, 2009.
 * @param p
 */
int random_int_geometric_from0(double p)
{
    if (p == 1)
        return 1;
    
    // generate random number
    double dum;
    do
        dum = random_double();
    while (dum == 0.0);
    
    // return random int
    return int(log(dum)/log(1-p));
}

/**
 * geometric random number generation (starting from 1)
 * Modified from W. Fletcher and Z. Yang, “INDELible: A flexible simulator of biological sequence evolution,” Mol. Biol. Evol., vol. 26, no. 8, pp. 1879–1888, 2009.
 * @param p
 */
int random_int_geometric(double p)
{
    return random_int_geometric_from0(p) + 1;
}

/**
 * negative binomial distribution
 * Modified from W. Fletcher and Z. Yang, “INDELible: A flexible simulator of biological sequence evolution,” Mol. Biol. Evol., vol. 26, no. 8, pp. 1879–1888, 2009.
 * @param r, q
 */
int random_int_nebin(int r, double q)
{
    int u = 0;
    while ( r-- )
        u += random_int_geometric_from0(1-q);
    return u + 1;
}

/**
 * Zipfian distribution
 * algorithm from DAWG (Cartwright, 2005).
 * Draw from Zipf distribution, with parameter a > 1.0
 * Devroye Luc (1986) Non-uniform random variate generation.
 * Springer-Verlag: Berlin. p551
 * @param a, m
 */
int random_int_zipf(double a, int m)
{
    double x;
    for (int i = 0; i < 1000; i++)
    {
        double b = pow(2.0, a-1.0);
        double t;
        do {
         x = floor(pow(random_double(), -1.0/(a-1.0)));
         t = pow(1.0+1.0/x, a-1.0);
        } while( random_double()*x*(t-1.0)*b > t*(b-1.0));
        
        // make sure x is not greater than the maximum value if m is set
        if (m == -1 || x <= m)
            break;
    }
    
    if (m != -1 && x > m)
        outError("Cannot generate a random length which is not greater than "+convertIntToString(m)+" based on the Zipfian distribution with the parameter a = "+convertDoubleToString(a)+". Please change either the parameter or the maximum value and try again!");
    return (int)x;
}

/**
 * Lavalette distribution
 * Modified from W. Fletcher and Z. Yang, “INDELible: A flexible simulator of biological sequence evolution,” Mol. Biol. Evol., vol. 26, no. 8, pp. 1879–1888, 2009.
 * @param a, m
 */
int random_int_lav(double a, int m)
{
    // initialize totald_vec
    double totald = 0;
    vector<double> totald_vec;
    for (double u = 1; u < m+1; u++)
    {
        totald += pow(m*u/(m-u+1),-1*a);
        totald_vec.push_back(totald);
    }
    
    // normalize totald_vec
    for(int i = 0; i < totald_vec.size(); i++)
    {
        totald_vec[i] /= totald;
    }
    
    // init random int
    double random_num = random_double();
    for(int i = 0; i < totald_vec.size(); i++)
        if (random_num < totald_vec.at(i))
            return i+1;

    return totald_vec.size();
}

/**
 * Parse indel-size distribution
 */
IndelDistribution parseIndelDis(string input, string event_name)
{
    // by default, using the user_defined distribution
    IndelDistribution indel_dis = IndelDistribution(USER_DEFINED, -1, -1, input);
    
    // detect the seperator
    char delimiter = ',';
    if (input.find('/') != std::string::npos)
        delimiter = '/';
    size_t pos;
    
    // parse negative binomial distribution
    if (input.rfind("nb{", 0) == 0 || input.rfind("NB{", 0) == 0)
    {
        // remove "nb{"
        input.erase(0, 3);
        
        // validate the parameters
        if (input[input.length()-1]!='}')
            throw "Use NB{<mean>[/<variance>]}";
        
        // remove "}"
        input = input.substr(0, input.length()-1);
        
        // determine the position of the delimiter (if any)
        pos = input.find(delimiter);
        
        // default value of r is 1
        int r = 1;
        // get mean and variance
        double mean, variance;
        // if users supply only a mean
        if (pos == std::string::npos)
        {
            // get mean
            mean = convert_double(input.c_str());
            
            // validate mean
            if (mean < 1)
                throw "<mean> must not less than 1";
        }
        // if users supply a mean and a variance
        else
        {
            // get mean
            mean = convert_double(input.substr(0, pos).c_str());
            
            // validate mean
            if (mean < 1)
                throw "<mean> must not less than 1";

            // get variance
            input.erase(0, pos + 1);
            double variance = convert_double(input.c_str());
            if (variance <= mean - 1)
                throw "<variance> must be greater than mean - 1";
            
            // compute r = (mean-1)^2/(variance - mean + 1)
            r = round((mean-1)*(mean-1)/(variance - mean + 1));
        }
        
        // compute q = (mean - 1)/(r + mean - 1)
        double q = (mean - 1)/(r + mean - 1);
        
        // show infor
        Params::getInstance().delay_msgs += "INFO: " + event_name + "-size is generated from Negative Binomial distribution with a mean of " + convertDoubleToString(mean) + " and a variance of " +convertDoubleToString((mean - 1)/(1 - q)) + ". The variance may be different from the user-defined value (if any) due to the definition of each distribution.\n";
        
        // return indel_dis
        return IndelDistribution(NEG_BIN, r, q);
    }
    
    // parse zipf
    if (input.rfind("pow", 0) == 0 || input.rfind("POW", 0) == 0)
    {
        // remove "pow"
        input.erase(0, 3);
        
        // determine the position of the delimiter
        pos = input.find(delimiter);
        
        // validate the parameters
        if ((input[0]!='{')
            || (input[input.length()-1]!='}'))
            throw "Use POW{<double_a>[/<int_max>]}";
        
        // get a
        double a;
        if (pos!= std::string::npos)
            a = convert_double(input.substr(1, pos - 1).c_str());
        else
            a = convert_double(input.substr(1, input.length() - 2).c_str());
        if (a <= 1)
            throw "<double_a> must be greater than 1";
        
        // get max (if any)
        int max = -1;
        if (pos != std::string::npos)
        {
            input.erase(0, pos + 1);
            max = convert_int(input.substr(0, input.length()-1).c_str());
            if (max <= 0)
                throw "<int_max> must be a positive integer";
        }
        
        // show infor
        string msg = "INFO: " + event_name + "-size is generated from Zipfian distribution with parameter <a> of " + convertDoubleToString(a);
        if (max != -1)
            msg += " and <max> of " + convertDoubleToString(max);
        msg += ".\n";
        
        Params::getInstance().delay_msgs += msg;
        
        // return indel_dis
        return IndelDistribution(ZIPF, a, max);
    }
    
    // parse Lavalette distribution
    if (input.rfind("lav", 0) == 0 || input.rfind("LAV", 0) == 0)
    {
        // remove "lav"
        input.erase(0, 3);
        
        // determine the position of the delimiter
        pos = input.find(delimiter);
        
        // validate the parameters
        if ((input[0]!='{')
            || (input[input.length()-1]!='}')
            || (pos == std::string::npos))
            throw "Use LAV{<double_a>/<int_max>}";
        
        // get a
        double a = convert_double(input.substr(1, pos - 1).c_str());
        
        // get max
        input.erase(0, pos + 1);
        int max = convert_int(input.substr(0, input.length()-1).c_str());
        if (max <= 0)
            throw "<int_max> must be a positive integer";
        
        // show infor
        Params::getInstance().delay_msgs += "INFO: " + event_name + "-size is generated from Lavalette distribution with parameter <a> of " + convertDoubleToString(a) + " and <max> of " + convertDoubleToString(max) + ".\n";
        
        // return indel_dis
        return IndelDistribution(LAV, a, max);
    }
    
    // parse Geometric distribution
    if (input.rfind("geo{", 0) == 0 || input.rfind("GEO{", 0) == 0)
    {
        // remove "geo{"
        input.erase(0, 4);
        
        // validate the parameters
        if (input[input.length()-1]!='}')
            throw "Use GEO{<mean>}";
        
        // remove "}"
        input = input.substr(0, input.length()-1);
            
        // determine the position of the delimiter (if any)
        pos = input.find(delimiter);
        
        // show warning if users specify a variance
        if (pos != std::string::npos)
        {
            Params::getInstance().delay_msgs += "In Geometric distribution, the variance could be computed from the mean. Thus, we ignore the user-specified variance.\n";
            
            // remove variance from the input
            input = input.substr(0, pos);
        }
        
        // get mean
        double mean = convert_double(input.c_str());
        if (mean < 1)
            throw "<mean> must not less than 1";
        // convert mean into p (of the distribution)
        double p = 1/mean;
        
        // show infor
        Params::getInstance().delay_msgs += "INFO: " + event_name + "-size is generated from Geometric distribution with a mean of " + convertDoubleToString(mean) + " and a variance of " + convertDoubleToString((1-p)/(p*p)) + ". The variance may be different from the user-defined value (if any) due to the definition of each distribution. \n";
        
        // return indel_dis
        return IndelDistribution(GEO, p);
    }
    
    return indel_dis;
}

void random_resampling(int n, IntVector &sample, int *rstream) {
    sample.resize(n, 0);
    if (Params::getInstance().jackknife_prop == 0.0) {
        // boostrap resampling
        for (int i = 0; i < n; i++) {
            int j = random_int(n, rstream);
            sample[j]++;
        }
    } else {
        // jackknife resampling
        int total = floor((1.0 - Params::getInstance().jackknife_prop)*n);
        if (total <= 0)
            outError("Jackknife sample size is zero");
        // make sure jackknife samples have exacly the same size
        for (int num = 0; num < total; ) {
            for (int i = 0; i < n; i++) if (!sample[i]) {
                if (random_double(rstream) < Params::getInstance().jackknife_prop)
                    continue;
                sample[i] = 1;
                num++;
                if (num >= total)
                    break;
            }
        }
    }
}


/* Following part is taken from ModelTest software */
#define	BIGX            20.0                                 /* max value to represent exp (x) */
#define	LOG_SQRT_PI     0.5723649429247000870717135          /* log (sqrt (pi)) */
#define	I_SQRT_PI       0.5641895835477562869480795          /* 1 / sqrt (pi) */
#define	Z_MAX           6.0                                  /* maximum meaningful z value */
#define	ex(x)           (((x) < -BIGX) ? 0.0 : exp (x))

/************** Normalz: probability of normal z value *********************/

/*
ALGORITHM:	Adapted from a polynomial approximation in:
                        Ibbetson D, Algorithm 209
                        Collected Algorithms of the CACM 1963 p. 616
                Note:
                        This routine has six digit accuracy, so it is only useful for absolute
                        z values < 6.  For z values >= to 6.0, Normalz() returns 0.0.
 */

double Normalz(double z) /*VAR returns cumulative probability from -oo to z VAR normal z value */ {
    double y, x, w;

    if (z == 0.0)
        x = 0.0;
    else {
        y = 0.5 * fabs(z);
        if (y >= (Z_MAX * 0.5))
            x = 1.0;
        else if (y < 1.0) {
            w = y*y;
            x = ((((((((0.000124818987 * w
                    - 0.001075204047) * w + 0.005198775019) * w
                    - 0.019198292004) * w + 0.059054035642) * w
                    - 0.151968751364) * w + 0.319152932694) * w
                    - 0.531923007300) * w + 0.797884560593) * y * 2.0;
        } else {
            y -= 2.0;
            x = (((((((((((((-0.000045255659 * y
                    + 0.000152529290) * y - 0.000019538132) * y
                    - 0.000676904986) * y + 0.001390604284) * y
                    - 0.000794620820) * y - 0.002034254874) * y
                    + 0.006549791214) * y - 0.010557625006) * y
                    + 0.011630447319) * y - 0.009279453341) * y
                    + 0.005353579108) * y - 0.002141268741) * y
                    + 0.000535310849) * y + 0.999936657524;
        }
    }
    return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}


/**************  ChiSquare: probability of chi square value *************/

/*ALGORITHM Compute probability of chi square value.
Adapted from: 	Hill, I. D. and Pike, M. C.  Algorithm 299.Collected Algorithms for the CACM 1967 p. 243
Updated for rounding errors based on remark inACM TOMS June 1985, page 185. Found in Perlman.lib*/

double computePValueChiSquare(double x, int df) /* x: obtained chi-square value,  df: degrees of freedom */ {
    double a, y, s;
    double e, c, z;
    int even; /* true if df is an even number */

    if (x <= 0.0 || df < 1)
        return (1.0);

    y = 1;

    a = 0.5 * x;
    even = (2 * (df / 2)) == df;
    if (df > 1)
        y = ex(-a);
    s = (even ? y : (2.0 * Normalz(-sqrt(x))));
    if (df > 2) {
        x = 0.5 * (df - 1.0);
        z = (even ? 1.0 : 0.5);
        if (a > BIGX) {
            e = (even ? 0.0 : LOG_SQRT_PI);
            c = log(a);
            while (z <= x) {
                e = log(z) + e;
                s += ex(c * z - a - e);
                z += 1.0;
            }
            return (s);
        } else {
            e = (even ? 1.0 : (I_SQRT_PI / sqrt(a)));
            c = 0.0;
            while (z <= x) {
                e = e * (a / z);
                c = c + e;
                z += 1.0;
            }
            return (c * y + s);
        }
    } else
        return (s);
}

void trimString(string &str) {
    str.erase(0, str.find_first_not_of(" \n\r\t"));
    str.erase(str.find_last_not_of(" \n\r\t")+1);
}



Params& Params::getInstance() {
    static Params instance;
    return instance;
}

void Params::setDefault() {
    tree_gen = NONE;
    user_file = nullptr;
    constraint_tree_file = nullptr;
    opt_gammai = true;
    opt_gammai_fast = false;
    opt_gammai_keep_bran = false;
    testAlphaEpsAdaptive = false;
    randomAlpha = false;
    testAlphaEps = 0.1;
    exh_ai = false;
    alpha_invar_file = nullptr;
    out_prefix = nullptr;
    out_file = nullptr;
    sub_size = 4;
    pd_proportion = 0.0;
    min_proportion = 0.0;
    step_proportion = 0.01;
    min_size = 0;
    step_size = 1;
    find_all = false;
    run_mode = RunMode::DETECTED;
    detected_mode = RunMode::DETECTED;
    param_file = nullptr;
    initial_file = nullptr;
    initial_area_file = nullptr;
    pdtaxa_file = nullptr;
    areas_boundary_file = nullptr;
    boundary_modifier = 1.0;
    dist_file = nullptr;
    compute_obs_dist = false;
    compute_jc_dist = true;
    experimental = true;
    compute_ml_dist = true;
    compute_ml_tree = true;
    compute_ml_tree_only = false;
    budget_file = nullptr;
    overlap = 0;
    is_rooted = false;
    root_move_dist = 2;
    root_find = false;
    root_test = false;
    sample_size = -1;
    repeated_time = 1;
    //nr_output = 10000;
    nr_output = 0;
    //smode = EXHAUSTIVE;
    intype = IN_OTHER;
    budget = -1;
    min_budget = -1;
    step_budget = 1;
    root = nullptr;
    num_splits = 0;
    min_len = 0.001;
    mean_len = 0.1;
    max_len = 0.999;
    num_zero_len = 0;
    pd_limit = 100;
    calc_pdgain = false;
    multi_tree = false;
    second_tree = nullptr;
    support_tag = nullptr;
    site_concordance = 0;
    ancestral_site_concordance = 0;
    site_concordance_partition = false;
    print_cf_quartets = false;
    print_df1_trees = false;
    internode_certainty = 0;
    tree_weight_file = nullptr;
    consensus_type = CT_NONE;
    find_pd_min = false;
    branch_cluster = 0;
    taxa_order_file = nullptr;
    endemic_pd = false;
    exclusive_pd = false;
    complement_area = nullptr;
    scaling_factor = -1;
    numeric_precision = -1;
    binary_programming = false;
    quad_programming = false;
    test_input = TEST_NONE;
    tree_burnin = 0;
    tree_max_count = 1000000;
    split_threshold = 0.0;
    split_threshold_str = nullptr;
    split_weight_threshold = -1000;
    collapse_zero_branch = false;
    split_weight_summary = SW_SUM;
    gurobi_format = true;
    gurobi_threads = 1;
    num_bootstrap_samples = 0;
    bootstrap_spec = nullptr;
    transfer_bootstrap = 0;

    aln_file = nullptr;
    phylip_sequential_format = false;
    symtest = SYMTEST_NONE;
    symtest_only = false;
    symtest_remove = 0;
    symtest_keep_zero = false;
    symtest_type = 0;
    symtest_pcutoff = 0.05;
    symtest_stat = false;
    symtest_shuffle = 1;
    //treeset_file = nullptr;
    topotest_replicates = 0;
    topotest_optimize_model = false;
    do_weighted_test = false;
    do_au_test = false;
    siteLL_file = nullptr; //added by MA
    partition_file = nullptr;
    partition_type = BRLEN_OPTIMIZE;
    partfinder_rcluster = 10; // change the default from 100 to 10
    partfinder_rcluster_max = 0;
    partition_merge = MERGE_NONE;
    merge_models = "1";
    merge_rates = "1";
    partfinder_log_rate = true;
    
    sequence_type = nullptr;
    aln_output = nullptr;
    aln_site_list = nullptr;
    aln_output_format = IN_PHYLIP;
    output_format = FORMAT_NORMAL;
    newick_extended_format = false;
    gap_masked_aln = nullptr;
    concatenate_aln = nullptr;
    aln_nogaps = false;
    aln_no_const_sites = false;
    print_aln_info = false;
//    parsimony = false;
//    parsimony_tree = false;
    tree_spr = false;
    nexus_output = false;
    k_representative = 4;
    loglh_epsilon = 0.001;
    numSmoothTree = 1;
    nni5 = true;
    nni5_num_eval = 1;
    brlen_num_traversal = 1;
    leastSquareBranch = false;
    pars_branch_length = false;
    bayes_branch_length = false;
    manuel_analytic_approx = false;
    leastSquareNNI = false;
    ls_var_type = OLS;
    maxCandidates = 20;
    popSize = 5;
    p_delete = -1;
    min_iterations = -1;
    max_iterations = 1000;
    num_param_iterations = 100;
    stop_condition = SC_UNSUCCESS_ITERATION;
    stop_confidence = 0.95;
    num_runs = 1;
    model_name = "";
    contain_nonrev = false;
    model_name_init = nullptr;
    model_opt_steps = 10;
    model_set = "ALL";
    model_extra_set = nullptr;
    model_subset = nullptr;
    state_freq_set = nullptr;
    ratehet_set = "AUTO";
    score_diff_thres = 10.0;
    model_def_file = nullptr;
    modelomatic = false;
    model_test_again = false;
    model_test_and_tree = 0;
    model_test_separate_rate = false;
    optimize_mixmodel_weight = false;
    optimize_mixmodel_freq = false;
    optimize_rate_matrix = false;
    store_trans_matrix = false;
    parallel_over_sites = false;
    order_by_threads = false;
    //freq_type = FREQ_EMPIRICAL;
    freq_type = FREQ_UNKNOWN;
    keep_zero_freq = true;
    min_state_freq = MIN_FREQUENCY;
    min_rate_cats = 2;
    num_rate_cats = 4;
    max_rate_cats = 10;
    min_mix_cats = 1;
    max_mix_cats = 10;
    start_subst = "GTR+FO";
    opt_rhas_again = false;
    opt_qmix_criteria = 2; // 1 : likelihood-ratio test; 2 : information criteria, like AIC, BIC
    opt_qmix_pthres = 0.05;
    check_combin_q_mat = true;
    gamma_shape = -1.0;
    min_gamma_shape = MIN_GAMMA_SHAPE;
    gamma_median = false;
    p_invar_sites = -1.0;
    optimize_model_rate_joint = false;
    optimize_by_newton = true;
    optimize_alg_freerate = "2-BFGS,EM";
    optimize_alg_mixlen = "EM";
    optimize_alg_gammai = "EM";
    optimize_alg_treeweight = "EM";
    optimize_from_given_params = false;
    optimize_alg_qmix = "BFGS";
    estimate_init_freq = 0;

    // defaults for new options -JD
    optimize_linked_gtr = false;
    gtr20_model = "POISSON";
    guess_multiplier = 0.75; // change from 0.5
    // rates_file = false;
    reset_method = "random"; // change from const

    optimize_params_use_hmm = false;
    optimize_params_use_hmm_sm = false;
    optimize_params_use_hmm_gm = false;
    optimize_params_use_hmm_tm = false;
    HMM_no_avg_brlen = false;
    HMM_min_stran = 0.0;
    treemix_optimize_methods = "mast"; // default is MAST

    fixed_branch_length = BRLEN_OPTIMIZE;
    min_branch_length = 0.0; // this is now adjusted later based on alignment length
    // TODO DS: This seems inappropriate for PoMo.  It is handled in
    // phyloanalysis::2908.
    max_branch_length = 10.0; // Nov 22 2016: reduce from 100 to 10!
    iqp_assess_quartet = IQP_DISTANCE;
    iqp = false;
    write_intermediate_trees = 0;
//    avoid_duplicated_trees = false;
    writeDistImdTrees = false;
    rf_dist_mode = 0;
    rf_same_pair = false;
    normalize_tree_dist = false;
    mvh_site_rate = false;
    rate_mh_type = true;
    discard_saturated_site = false;
    mean_rate = 1.0;
    aLRT_threshold = 101;
    aLRT_replicates = 0;
    aLRT_test = false;
    aBayes_test = false;
    localbp_replicates = 0;
#ifdef __AVX512KNL
    SSE = LK_AVX512;
#else
    SSE = LK_AVX_FMA;
#endif
    lk_safe_scaling = false;
    numseq_safe_scaling = 2000;
    kernel_nonrev = false;
    print_site_lh = WSL_NONE;
    print_partition_lh = false;
    print_marginal_prob = false;
    print_site_prob = WSL_NONE;
    print_site_state_freq = WSF_NONE;
    print_site_rate = 0;
    print_trees_site_posterior = 0;
    print_ancestral_sequence = AST_NONE;
    min_ancestral_prob = 0.0;
    print_tree_lh = false;
    lambda = 1;
    speed_conf = 1.0;
    whtest_simulations = 1000;
    mcat_type = MCAT_LOG + MCAT_PATTERN;
    rate_file = nullptr;
    ngs_file = nullptr;
    ngs_mapped_reads = nullptr;
    ngs_ignore_gaps = true;
    do_pars_multistate = false;
    gene_pvalue_file = nullptr;
    gene_scale_factor = -1;
    gene_pvalue_loga = false;
    second_align = nullptr;
    ncbi_taxid = 0;
    ncbi_taxon_level = nullptr;
    ncbi_names_file = nullptr;
    ncbi_ignore_level = nullptr;
    eco_dag_file  = nullptr;
    eco_type = nullptr;
    eco_detail_file = nullptr;
    k_percent = 0;
    diet_min = 0;
    diet_max = 0;
    diet_step = 0;
    eco_weighted = false;
    eco_run = 0;

    upper_bound = false;
    upper_bound_NNI = false;
    upper_bound_frac = 0.0;

    gbo_replicates = 0;
    ufboot_epsilon = 0.5;
    check_gbo_sample_size = 0;
    use_rell_method = true;
    use_elw_method = false;
    use_weighted_bootstrap = false;
    use_max_tree_per_bootstrap = true;
    max_candidate_trees = 0;
    distinct_trees = false;
    online_bootstrap = true;
    min_correlation = 0.99;
    step_iterations = 100;
//    store_candidate_trees = false;
    print_ufboot_trees = 0;
    jackknife_prop = 0.0;
    robust_phy_keep = 1.0;
    robust_median = false;
    //const double INF_NNI_CUTOFF = -1000000.0;
    nni_cutoff = -1000000.0;
    estimate_nni_cutoff = false;
    nni_sort = false;
    //nni_opt_5branches = false;
    testNNI = false;
    approximate_nni = false;
    do_compression = false;

    new_heuristic = true;
    iteration_multiple = 1;
    initPS = 0.5;
#ifdef USING_PLL
    pll = true;
#else
    pll = false;
#endif
    modelEps = 0.01;
    fundiEps = 0.000001;
    modelfinder_eps = 0.1;
    treemix_eps = 0.001;
    treemixhmm_eps = 0.01;
    parbran = false;
    binary_aln_file = nullptr;
    maxtime = 1000000;
    reinsert_par = false;
    bestStart = true;
    snni = true; // turn on sNNI default now
//    autostop = true; // turn on auto stopping rule by default now
    unsuccess_iteration = 100;
    speednni = true; // turn on reduced hill-climbing NNI by default now
    numInitTrees = 100;
    fixStableSplits = false;
    stableSplitThreshold = 0.9;
    five_plus_five = false;
    memCheck = false;
    tabu = false;
    adaptPertubation = false;
    numSupportTrees = 20;
//    sprDist = 20;
    sprDist = 6;
    sankoff_cost_file = nullptr;
    numNNITrees = 20;
    avh_test = 0;
    bootlh_test = 0;
    bootlh_partitions = nullptr;
    site_freq_file = nullptr;
    tree_freq_file = nullptr;
    num_threads = 1;
    num_threads_max = 10000;
    openmp_by_model = false;
    model_test_criterion = MTC_BIC;
//    model_test_stop_rule = MTC_ALL;
    model_test_sample_size = 0;
    root_state = nullptr;
    print_bootaln = false;
    print_boot_site_freq = false;
    print_subaln = false;
    print_partition_info = false;
    print_conaln = false;
    count_trees = false;
    pomo = false;
    pomo_random_sampling = false;
    // pomo_counts_file_flag = false;
    pomo_pop_size = 9;
    print_branch_lengths = false;
    lh_mem_save = LM_PER_NODE; // auto detect
    buffer_mem_save = false;
    start_tree = STT_PLL_PARSIMONY;
    start_tree_subtype_name = StartTree::Factory::getNameOfDefaultTreeBuilder();

    modelfinder_ml_tree = true;
    final_model_opt = true;
    print_splits_file = false;
    print_splits_nex_file = true;
    ignore_identical_seqs = true;
    write_init_tree = false;
    write_candidate_trees = false;
    write_branches = false;
    freq_const_patterns = nullptr;
    no_rescale_gamma_invar = false;
    compute_seq_identity_along_tree = false;
    compute_seq_composition = true;
    lmap_num_quartets = -1;
    lmap_cluster_file = nullptr;
    print_lmap_quartet_lh = false;
    num_mixlen = 1;
    link_alpha = false;
    link_model = false;
    model_joint = "";
    ignore_checkpoint = false;
    checkpoint_dump_interval = 60;
    force_unfinished = false;
    print_all_checkpoints = false;
    suppress_output_flags = 0;
    ufboot2corr = false;
    u2c_nni5 = false;
    date_with_outgroup = true;
    date_debug = false;
    date_replicates = 0;
    clock_stddev = -1.0;
    date_outlier = -1.0;
    dating_mf = false;
    mcmc_clock = CORRELATED;
    mcmc_bds = "1 1 0.5";
    mcmc_iter = "20000, 100, 20000";

    // added by TD
    use_nn_model = false;
    nn_path_model = "resnet_modelfinder.onnx";
    nn_path_rates = "lanfear_alpha_lstm.onnx";

    // ------------ Terrace variables ------------
    terrace_check = false;
    terrace_analysis = false;
    print_terrace_trees = false;
    print_induced_trees = false;
    pr_ab_matrix = nullptr;
    print_pr_ab_matrix = false;
    print_m_overlap = false;
    terrace_query_set = nullptr;
    terrace_stop_intermediate_num = -1;
    terrace_stop_terrace_trees_num = -1;
    terrace_stop_time = -1;
    terrace_non_stop = false;
    terrace_print_lim = 0;
    terrace_remove_m_leaves = 0;
    matrix_order = false;
    gen_all_NNI = false;
    
    remove_empty_seq = true;
    terrace_aware = true;
#ifdef IQTREE_TERRAPHAST
    terrace_analysis_tphast = false;
#else
    terrace_analysis_tphast = false;
#endif
    
    // --------------------------------------------
    
    matrix_exp_technique = MET_EIGEN3LIB_DECOMPOSITION;

    if (nni5) {
        nni_type = NNI5;
    } else {
        nni_type = NNI1;
    }

    struct timeval tv;
    struct timezone tz;
    // initialize random seed based on current time
    gettimeofday(&tv, &tz);
    //ran_seed = (unsigned) (tv.tv_sec+tv.tv_usec);
    ran_seed = (tv.tv_usec);
    subsampling_seed = ran_seed;
    subsampling = 0;
    
    suppress_list_of_sequences = false;
    suppress_zero_distance_warnings = false;
    suppress_duplicate_sequence_warnings = false;
    
    original_params = "";
    alisim_active = false;
    multi_rstreams_used = false;
    alisim_inference_mode = false;
    alisim_no_copy_gaps = false;
    alisim_sequence_length = 1000;
    alisim_dataset_num = 1;
    root_ref_seq_aln = "";
    root_ref_seq_name = "";
    alisim_max_rate_categories_for_applying_caching = 100;
    alisim_num_states_morph = 0;
    alisim_num_taxa_uniform_start = -1;
    alisim_num_taxa_uniform_end = -1;
    alisim_length_ratio = 2;
    birth_rate = 0.8;
    death_rate = 0.2;
    alisim_fundi_proportion = 0.0;
    fundi_init_proportion = 0.5;
    fundi_init_branch_length = 0.0;
    alisim_distribution_definitions = nullptr;
    alisim_skip_checking_memory = false;
    alisim_write_internal_sequences = false;
    alisim_only_unroot_tree = false;
    branch_distribution = nullptr;
    alisim_insertion_ratio = 0;
    alisim_deletion_ratio = 0;
    alisim_insertion_distribution = IndelDistribution(ZIPF,1.7,100);
    alisim_deletion_distribution = IndelDistribution(ZIPF,1.7,100);
    alisim_mean_deletion_size = -1;
    alisim_simulation_thresh = 0.001;
    delay_msgs = "";
    alisim_no_export_sequence_wo_gaps = false;
    alisim_mixture_at_sub_level = false;
    alisim_branch_scale = 1.0;
    alisim_rate_heterogeneity = POSTERIOR_MEAN;
    alisim_stationarity_heterogeneity = POSTERIOR_MEAN;
    alisim_single_output = false;
    keep_seq_order = false;
    mem_limit_factor = 0;
    delete_output = false;
    indel_rate_variation = false;
    tmp_data_filename = "tmp_data";
    rebuild_indel_history_param = 1.0/3;
    alisim_openmp_alg = IM;
    no_merge = false;
    alignment_id = 0;
    inference_alg = ALG_IQ_TREE;
    in_aln_format_str = "AUTO";
    shallow_tree_search = false;
    tree_search_type_str = "NORMAL";
    allow_replace_input_tree = false;
    tree_format_str = "BIN";
    make_consistent = false;
    include_pre_mutations = false;
    mutation_file = "";
    site_starting_index = 0;
    
    // ----------- SPRTA ----------
    compute_SPRTA = false;
    SPRTA_zero_branches = false;
    out_alter_spr = false;
    intree_str = "";
}

int countPhysicalCPUCores() {
    #ifdef _OPENMP
    return omp_get_num_procs();
    #else
    return std::thread::hardware_concurrency();
    #endif
    /*
    uint32_t registers[4];
    unsigned logicalcpucount;
    unsigned physicalcpucount;
#if defined(_WIN32) || defined(WIN32)
    SYSTEM_INFO systeminfo;
    GetSystemInfo( &systeminfo );
    logicalcpucount = systeminfo.dwNumberOfProcessors;
#else
    logicalcpucount = sysconf( _SC_NPROCESSORS_ONLN );
#endif
    if (logicalcpucount < 1) logicalcpucount = 1;
    return logicalcpucount;
    
    if (logicalcpucount % 2 != 0)
        return logicalcpucount;
    __asm__ __volatile__ ("cpuid " :
                          "=a" (registers[0]),
                          "=b" (registers[1]),
                          "=c" (registers[2]),
                          "=d" (registers[3])
                          : "a" (1), "c" (0));

    unsigned CPUFeatureSet = registers[3];
    bool hyperthreading = CPUFeatureSet & (1 << 28);    
    if (hyperthreading){
        physicalcpucount = logicalcpucount / 2;
    } else {
        physicalcpucount = logicalcpucount;
    }
    if (physicalcpucount < 1) physicalcpucount = 1;
    return physicalcpucount;
     */
}

// stacktrace.h (c) 2008, Timo Bingmann from http://idlebox.net/
// published under the WTFPL v2.0

/** Print a demangled stack backtrace of the caller function to FILE* out. */

#if  !defined(Backtrace_FOUND)

// donothing for WIN32
void print_stacktrace(ostream &out, unsigned int max_frames) {}

#else

void print_stacktrace(ostream &out, unsigned int max_frames)
{
#ifdef _OPENMP
#pragma omp master
{
#endif
    out << "STACK TRACE FOR DEBUGGING:" << endl;

    // storage array for stack trace address data
    void* addrlist[max_frames+1];

    // retrieve current stack addresses
    int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));

//    if (addrlen == 0) {
//        out << "  <empty, possibly corrupt>" << endl;
//        return;
//    }

    // resolve addresses into strings containing "filename(function+address)",
    // this array must be free()-ed
    char** symbollist = backtrace_symbols(addrlist, addrlen);

    // allocate string which will be filled with the demangled function name
    size_t funcnamesize = 256;
    char* funcname = (char*)malloc(funcnamesize);

    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for (int i = 1; i < addrlen; i++)
    {
	char *begin_name = 0, *begin_offset = 0;

	// find parentheses and +address offset surrounding the mangled name:
#ifdef __clang__
      // OSX style stack trace
      for ( char *p = symbollist[i]; *p; ++p )
      {
         if (( *p == '_' ) && ( *(p-1) == ' ' ))
            begin_name = p-1;
         else if ( *p == '+' )
            begin_offset = p-1;
      }

      if ( begin_name && begin_offset && ( begin_name < begin_offset ))
      {
         *begin_name++ = '\0';
         *begin_offset++ = '\0';

         // mangled name is now in [begin_name, begin_offset) and caller
         // offset in [begin_offset, end_offset). now apply
         // __cxa_demangle():
         int status;
         char* ret = abi::__cxa_demangle( begin_name, &funcname[0],
                                          &funcnamesize, &status );
         if ( status == 0 )
         {
            funcname = ret; // use possibly realloc()-ed string
//            out << "  " << symbollist[i] << " : " << funcname << "+"<< begin_offset << endl;
            out << i << "   "  << funcname << endl;
         } else {
            // demangling failed. Output function name as a C function with
            // no arguments.
//             out << "  " << symbollist[i] << " : " << begin_name << "()+"<< begin_offset << endl;
            out << i << "   " << begin_name << "()" << endl;
         }

#else // !DARWIN - but is posix
         // ./module(function+0x15c) [0x8048a6d]
    char *end_offset = 0;
	for (char *p = symbollist[i]; *p; ++p)
	{
	    if (*p == '(')
		begin_name = p;
	    else if (*p == '+')
		begin_offset = p;
	    else if (*p == ')' && begin_offset) {
		end_offset = p;
		break;
	    }
	}

	if (begin_name && begin_offset && end_offset
	    && begin_name < begin_offset)
	{
	    *begin_name++ = '\0';
	    *begin_offset++ = '\0';
	    *end_offset = '\0';

	    // mangled name is now in [begin_name, begin_offset) and caller
	    // offset in [begin_offset, end_offset). now apply
	    // __cxa_demangle():

	    int status;
	    char* ret = abi::__cxa_demangle(begin_name,
					    funcname, &funcnamesize, &status);
	    if (status == 0) {
            funcname = ret; // use possibly realloc()-ed string
//            out << "  " << symbollist[i] << " : " << funcname << "+"<< begin_offset << endl;
            out << i << "   " << funcname << endl;
	    }
	    else {
            // demangling failed. Output function name as a C function with
            // no arguments.
//            out << "  " << symbollist[i] << " : " << begin_name << "()+"<< begin_offset << endl;
            out << i << "   " << begin_name << "()" << endl;
	    }
#endif
	}
	else
	{
	    // couldn't parse the line? print the whole line.
//	    out << i << ". " << symbollist[i] << endl;
	}
    }

    free(funcname);
    free(symbollist);
#ifdef _OPENMP
}
#endif

}

#endif // Backtrace_FOUND

bool memcmpcpy(void * destination, const void * source, size_t num) {
    bool diff = (memcmp(destination, source, num) != 0);
    memcpy(destination, source, num);
    return diff;
}

// Pairing function: see https://en.wikipedia.org/wiki/Pairing_function
int pairInteger(int int1, int int2) {
    if (int1 <= int2) {
        return ((int1 + int2)*(int1 + int2 + 1)/2 + int2);
    } else {
        return ((int1 + int2)*(int1 + int2 + 1)/2 + int1);
    }
}

/*
 * Given a model name, look in it for "+F..." and 
 * determine the StateFreqType. Returns FREQ_EMPIRICAL if
 * unable to find a good +F... specifier
 */
StateFreqType parseStateFreqFromPlusF(string model_name) {
//    StateFreqType freq_type = FREQ_EMPIRICAL;

    // BQM 2017-05-02: change back to FREQ_UNKNOWN to resemble old behavior
    StateFreqType freq_type = FREQ_UNKNOWN;
    size_t plusFPos;
    if (model_name.find("+F1X4") != string::npos)
        freq_type = FREQ_CODON_1x4;
    else if (model_name.find("+F3X4C") != string::npos)
        freq_type = FREQ_CODON_3x4C;
    else if (model_name.find("+F3X4") != string::npos)
        freq_type = FREQ_CODON_3x4;
    else if (model_name.find("+FQ") != string::npos)
        freq_type = FREQ_EQUAL;
    else if (model_name.find("+FO") != string::npos)
        freq_type = FREQ_ESTIMATE;
    else if (model_name.find("+FU") != string::npos)
        freq_type = FREQ_USER_DEFINED;
    else if (model_name.find("+FRY") != string::npos)
        freq_type = FREQ_DNA_RY;
    else if (model_name.find("+FWS") != string::npos)
        freq_type = FREQ_DNA_WS;
    else if (model_name.find("+FMK") != string::npos)
        freq_type = FREQ_DNA_MK;
    else if ((plusFPos = model_name.find("+F")) != string::npos) {
        freq_type = FREQ_EMPIRICAL;
        // Now look for +F#### where #s are digits
        if (model_name.length() > plusFPos+2 && isdigit(model_name[plusFPos+2]))
        try {
            // throws if string is not 4 digits
            freq_type = parseStateFreqDigits(model_name.substr(plusFPos+2,4));
        } catch (const char *str) {
            // +F exists, but can't parse it as anything else
            outError(str);
        }
    }
    return(freq_type);
}

/*
 * Given a string of 4 digits, return a StateFreqType according to
 * equality constraints expressed by those digits.
 * E.g. "1233" constrains pi_G=pi_T (ACGT order, 3rd and 4th equal)
 * which results in FREQ_DNA_2311. "5288" would give the same result.
 */

StateFreqType parseStateFreqDigits(string digits) {
    bool good = true;
    if (digits.length()!=4) {
        good = false;
    } else {
        // Convert digits to canonical form, first occuring digit becomes 1 etc.
        int digit_order[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
        int first_found = 0;
        for (int i=0; i<4; i++) {
            int digit = digits[i]-'0';
            if (digit<0 || digit>9) {
	good = false; // found a non-digit
	break;
            }
            if (digit_order[digit]==-1) {
	// haven't seen this digit before
	digit_order[digit]=++first_found;
            }
            // rewrite digit in canonical form
            digits[i] = '0'+digit_order[digit];
        }
        // e.g. if digits was "5288", digit_order will end up as {-1,-1,2,-1,-1,1,-1,-1,3,-1}
    }
    if (!good) throw "Use -f <c | o | u | q | ry | ws | mk | <digit><digit><digit><digit>>";
    StateFreqType freq_type = FREQ_UNKNOWN;
    // Now just exhaustively list all canonical digit possibilities
    if (digits.compare("1111")==0) {
        freq_type = FREQ_EQUAL;
    } else if (digits.compare("1112")==0) {
        freq_type = FREQ_DNA_1112;
    } else if (digits.compare("1121")==0) {
        freq_type = FREQ_DNA_1121;
    } else if (digits.compare("1211")==0) {
        freq_type = FREQ_DNA_1211;
    } else if (digits.compare("1222")==0) {
        freq_type = FREQ_DNA_2111;
    } else if (digits.compare("1122")==0) {
        freq_type = FREQ_DNA_1122;
    } else if (digits.compare("1212")==0) {
        freq_type = FREQ_DNA_1212;
    } else if (digits.compare("1221")==0) {
        freq_type = FREQ_DNA_1221;
    } else if (digits.compare("1123")==0) {
        freq_type = FREQ_DNA_1123;
    } else if (digits.compare("1213")==0) {
        freq_type = FREQ_DNA_1213;
    } else if (digits.compare("1231")==0) {
        freq_type = FREQ_DNA_1231;
    } else if (digits.compare("1223")==0) {
        freq_type = FREQ_DNA_2113;
    } else if (digits.compare("1232")==0) {
        freq_type = FREQ_DNA_2131;
    } else if (digits.compare("1233")==0) {
        freq_type = FREQ_DNA_2311;
    } else if (digits.compare("1234")==0) {
        freq_type = FREQ_ESTIMATE;
    } else
        throw ("Unrecognized canonical digits - Can't happen"); // paranoia is good.
    return freq_type;
}


/*
 * All params in range [0,1] 
 * returns true if base frequencies have changed as a result of this call
 */

bool freqsFromParams(double *freq_vec, double *params, StateFreqType freq_type) {

    // BQM 2017-05-02: Note that only freq for A, C, G are free parameters and stored
    // in params, whereas freq_T is not free and should be handled properly

    double pA, pC, pG, pT; // base freqs
    switch (freq_type) {
    case FREQ_EQUAL:
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
        return false;
    case FREQ_ESTIMATE:
    // Minh: in code review, please pay extra attention to ensure my treadment of FREQ_ESTIMATE is equivalent to old treatment.
    // BQM: DONE!
        pA=params[0];
        pC=params[1];
        pG=params[2];
        //pT=1-pA-pC-pG;
        pT=freq_vec[3];
        break;
    case FREQ_DNA_RY:
        pA = params[0]/2;
        pG = 0.5-pA;
        pC = params[1]/2;
        pT = 0.5-pC;
        break;
    case FREQ_DNA_WS:
        pA = params[0]/2;
        pT = 0.5-pA;
        pC = params[1]/2;
        pG = 0.5-pC;
        break;
    case FREQ_DNA_MK:
        pA = params[0]/2;
        pC = 0.5-pA;
        pG = params[1]/2;
        pT = 0.5-pG;
        break;
    case FREQ_DNA_1112:
        pA = pC = pG = params[0]/3;
        pT = 1-3*pA;
        break;
    case FREQ_DNA_1121:
        pA = pC = pT = params[0]/3;
        pG = 1-3*pA;
        break;
    case FREQ_DNA_1211:
        pA = pG = pT = params[0]/3;
        pC = 1-3*pA;
        break;
    case FREQ_DNA_2111:
        pC = pG = pT = params[0]/3;
        pA = 1-3*pC;
        break;
    case FREQ_DNA_1122:
        pA = params[0]/2;
        pC = pA;
        pG = 0.5-pA;
        pT = pG;
        break;
    case FREQ_DNA_1212:
        pA = params[0]/2;
        pG = pA;
        pC = 0.5-pA;
        pT = pC;
        break;
    case FREQ_DNA_1221:
        pA = params[0]/2;
        pT = pA;
        pC = 0.5-pA;
        pG = pC;
        break;
    case FREQ_DNA_1123:
        pA = params[0]/2;
        pC = pA;
        pG = params[1]*(1-2*pA);
        pT = 1-pG-2*pA;
        break;
    case FREQ_DNA_1213:
        pA = params[0]/2;
        pG = pA;
        pC = params[1]*(1-2*pA);
        pT = 1-pC-2*pA;
        break;
    case FREQ_DNA_1231:
        pA = params[0]/2;
        pT = pA;
        pC = params[1]*(1-2*pA);
        pG = 1-pC-2*pA;
        break;
    case FREQ_DNA_2113:
        pC = params[0]/2;
        pG = pC;
        pA = params[1]*(1-2*pC);
        pT = 1-pA-2*pC;
        break;
    case FREQ_DNA_2131:
        pC = params[0]/2;
        pT = pC;
        pA = params[1]*(1-2*pC);
        pG = 1-pA-2*pC;
        break;
    case FREQ_DNA_2311:
        pG = params[0]/2;
        pT = pG;
        pA = params[1]*(1-2*pG);
        pC = 1-pA-2*pG;
        break;
    default:
        throw("Unrecognized freq_type in freqsFromParams - can't happen");
    }

    // To MDW, 2017-05-02: please make sure that frequencies are positive!
    // Otherwise, numerical issues will occur.

    bool changed = freq_vec[0]!=pA || freq_vec[1]!=pC || freq_vec[2]!=pG || freq_vec[3]!=pT;
    if (changed) {
        freq_vec[0]=pA;
        freq_vec[1]=pC;
        freq_vec[2]=pG;
        freq_vec[3]=pT;
    }
    return(changed);
}

/*
 * For given freq_type, derives frequency parameters from freq_vec
 * All parameters are in range [0,1] (assuming freq_vec is valid)
 */

void paramsFromFreqs(double *params, double *freq_vec, StateFreqType freq_type) {
    double pA = freq_vec[0]; // These just improve code readability
    double pC = freq_vec[1];
    double pG = freq_vec[2];
//    double pT = freq_vec[3]; // pT is not used below
    switch (freq_type) {
    case FREQ_EQUAL:
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
        break; // freq_vec never changes
    case FREQ_ESTIMATE:
        params[0]=pA;
        params[1]=pC;
        params[2]=pG;
        break;
    case FREQ_DNA_RY:
        params[0]=2*pA;
        params[1]=2*pC;
        break;
    case FREQ_DNA_WS:
        params[0]=2*pA;
        params[1]=2*pC;
        break;
    case FREQ_DNA_MK:
        params[0]=2*pA;
        params[1]=2*pG;
        break;
    case FREQ_DNA_1112:
        params[0]=3*pA;
        break;
    case FREQ_DNA_1121:
        params[0]=3*pA;
        break;
    case FREQ_DNA_1211:
        params[0]=3*pA;
        break;
    case FREQ_DNA_2111:
        params[0]=3*pC;
        break;
    case FREQ_DNA_1122:
        params[0]=2*pA;
        break;
    case FREQ_DNA_1212:
        params[0]=2*pA;
        break;
    case FREQ_DNA_1221:
        params[0]=2*pA;
        break;
    case FREQ_DNA_1123:
        params[0]=2*pA;
        params[1]=pG/(1-params[0]);
        break;
    case FREQ_DNA_1213:
        params[0]=2*pA;
        params[1]=pC/(1-params[0]);
        break;
    case FREQ_DNA_1231:
        params[0]=2*pA;
        params[1]=pC/(1-params[0]);
        break;
    case FREQ_DNA_2113:
        params[0]=2*pC;
        params[1]=pA/(1-params[0]);
        break;
    case FREQ_DNA_2131:
        params[0]=2*pC;
        params[1]=pA/(1-params[0]);
        break;
    case FREQ_DNA_2311:
        params[0]=2*pG;
        params[1]=pA/(1-params[0]);
        break;
    default:
        throw("Unrecognized freq_type in paramsFromFreqs - can't happen");
    }
}

/* 
 * Given a DNA freq_type and a base frequency vector, alter the
 * base freq vector to conform with the constraints of freq_type
 */
void forceFreqsConform(double *base_freq, StateFreqType freq_type) {
    double pA = base_freq[0]; // These just improve code readability
    double pC = base_freq[1];
    double pG = base_freq[2];
    double pT = base_freq[3];
    double scale;
    switch (freq_type) {
    case FREQ_EQUAL:
        // this was already handled, thus not necessary to check here 
//        base_freq[0] = base_freq[1] = base_freq[2] = base_freq[3] = 0.25;
//        break;
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
    case FREQ_ESTIMATE:
        break; // any base_freq is legal
    case FREQ_DNA_RY:
        scale = 0.5/(pA+pG);
        base_freq[0] = pA*scale;
        base_freq[2] = pG*scale;
        scale = 0.5/(pC+pT);
        base_freq[1] = pC*scale;
        base_freq[3] = pT*scale;
        break;
    case FREQ_DNA_WS:
        scale = 0.5/(pA+pT);
        base_freq[0] = pA*scale;
        base_freq[3] = pT*scale;
        scale = 0.5/(pC+pG);
        base_freq[1] = pC*scale;
        base_freq[2] = pG*scale;
        break;
    case FREQ_DNA_MK:
        scale = 0.5/(pA+pC);
        base_freq[0] = pA*scale;
        base_freq[1] = pC*scale;
        scale = 0.5/(pG+pT);
        base_freq[2] = pG*scale;
        base_freq[3] = pT*scale;
        break;
    case FREQ_DNA_1112:
        base_freq[0]=base_freq[1]=base_freq[2]=(pA+pC+pG)/3;
        break;
    case FREQ_DNA_1121:
        base_freq[0]=base_freq[1]=base_freq[3]=(pA+pC+pT)/3;
        break;
    case FREQ_DNA_1211:
        base_freq[0]=base_freq[2]=base_freq[3]=(pA+pG+pT)/3;
        break;
    case FREQ_DNA_2111:
        base_freq[1]=base_freq[2]=base_freq[3]=(pC+pG+pT)/3;
        break;
    case FREQ_DNA_1122:
        base_freq[0]=base_freq[1]=(pA+pC)/2;
        base_freq[2]=base_freq[3]=(pG+pT)/2;
        break;
    case FREQ_DNA_1212:
        base_freq[0]=base_freq[2]=(pA+pG)/2;
        base_freq[1]=base_freq[3]=(pC+pT)/2;
        break;
    case FREQ_DNA_1221:
        base_freq[0]=base_freq[3]=(pA+pT)/2;
        base_freq[1]=base_freq[2]=(pC+pG)/2;
        break;
    case FREQ_DNA_1123:
        base_freq[0]=base_freq[1]=(pA+pC)/2;
        break;
    case FREQ_DNA_1213:
        base_freq[0]=base_freq[2]=(pA+pG)/2;
        break;
    case FREQ_DNA_1231:
        base_freq[0]=base_freq[3]=(pA+pT)/2;
        break;
    case FREQ_DNA_2113:
        base_freq[1]=base_freq[2]=(pC+pG)/2;
        break;
    case FREQ_DNA_2131:
        base_freq[1]=base_freq[3]=(pC+pT)/2;
        break;
    case FREQ_DNA_2311:
        base_freq[2]=base_freq[3]=(pG+pT)/2;
        break;
    default:
        throw("Unrecognized freq_type in forceFreqsConform - can't happen");
    }
    ASSERT(base_freq[0]>=0 && base_freq[1]>=0 && base_freq[2]>=0 && base_freq[3]>=0 && fabs(base_freq[0]+base_freq[1]+base_freq[2]+base_freq[3]-1)<1e-7);
}

/*
 * For given freq_type, how many parameters are needed to
 * determine frequenc vector?
 * Currently, this is for DNA StateFreqTypes only.
 */

int nFreqParams(StateFreqType freq_type) {
    switch (freq_type) {
    case FREQ_DNA_1112:
    case FREQ_DNA_1121:
    case FREQ_DNA_1211:
    case FREQ_DNA_2111:
    case FREQ_DNA_1122:
    case FREQ_DNA_1212:
    case FREQ_DNA_1221:
        return(1);
    case FREQ_DNA_RY:
    case FREQ_DNA_WS:
    case FREQ_DNA_MK:
    case FREQ_DNA_1123:
    case FREQ_DNA_1213:
    case FREQ_DNA_1231:
    case FREQ_DNA_2113:
    case FREQ_DNA_2131:
    case FREQ_DNA_2311:
        return(2);   
    default:
        return 0; // BQM: don't care about other cases
    }
}

/*
 * For freq_type, and given every base must have frequency >= min_freq, set upper
 * and lower bounds for parameters.
 */
 void setBoundsForFreqType(double *lower_bound, 
                           double *upper_bound, 
                           bool *bound_check, 
                           double min_freq, 
                           StateFreqType freq_type) {
    // Sanity check: if min_freq==0, lower_bound=0 and upper_bound=1 
    // (except FREQ_ESTIMATE, which follows legacy code way of doing things.)
    switch (freq_type) {
    case FREQ_EQUAL:
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
        break; // There are no frequency determining parameters
    case FREQ_DNA_1112:
    case FREQ_DNA_1121:
    case FREQ_DNA_1211:
    case FREQ_DNA_2111:
        // one frequency determining parameter
        lower_bound[0] = 3*min_freq;
        upper_bound[0] = 1-min_freq;
        bound_check[0] = true;
        break;
    case FREQ_DNA_1122:
    case FREQ_DNA_1212:
    case FREQ_DNA_1221:
        // one frequency determining parameter
        lower_bound[0] = 2*min_freq;
        upper_bound[0] = 1-2*min_freq;
        bound_check[0] = true;
        break;
    case FREQ_DNA_RY:
    case FREQ_DNA_WS:
    case FREQ_DNA_MK:
        // two frequency determining parameters
        lower_bound[0] = lower_bound[1] = 2*min_freq;
        upper_bound[0] = upper_bound[1] = 1-2*min_freq;
        bound_check[0] = bound_check[1] = true;
	break;
    case FREQ_DNA_1123:
    case FREQ_DNA_1213:
    case FREQ_DNA_1231:
    case FREQ_DNA_2113:
    case FREQ_DNA_2131:
    case FREQ_DNA_2311:
        // two frequency determining parameters
        lower_bound[0] = 2*min_freq;
        upper_bound[0] = 1-2*min_freq;
	lower_bound[1] = min_freq/(1-2*min_freq);
        upper_bound[1] = (1-3*min_freq)/(1-2*min_freq);
        bound_check[0] = bound_check[1] = true;
	break;
        /* NOTE:
	 * upper_bound[1] and lower_bound[1] are not perfect. Some in-bounds parameters
         * will give base freqs for '2' or '3' base below minimum. This is
         * the best that can be done without passing min_freq to freqsFromParams
         */
    case FREQ_ESTIMATE:
        lower_bound[0] = lower_bound[1] = lower_bound[2] = min_freq;
        upper_bound[0] = upper_bound[1] = upper_bound[2] = 1;
        bound_check[0] = bound_check[1] = bound_check[2] = false;
        break;
    default:
        throw("Unrecognized freq_type in setBoundsForFreqType - can't happen");
    }
}
 
double binomial_coefficient_log(unsigned int N, unsigned int n) {
  static DoubleVector logv;
  if (logv.size() <= 0) {
    logv.push_back(0.0);
    logv.push_back(0.0);
  }
  if (n < N-n)
    n = N-n;
  if (n==0)
    return 0.0;
  if (N >= logv.size()) {
    for (unsigned int i = logv.size(); i <= N; i++)
      logv.push_back(log((double) i));
  }
  double binom_log = 0.0;
  for (unsigned int i = n+1; i <= N; i++)
    binom_log += logv[i] - logv[i-n];
  return binom_log;
}

double binomial_dist(unsigned int k, unsigned int N, double p) {
  double binom_log = binomial_coefficient_log(N, k);
  double res_log = binom_log + log(p)*k + log(1-p)*(N-k);
  return exp(res_log);
}

double hypergeometric_dist(unsigned int k, unsigned int n, unsigned int K, unsigned int N) {
  if (n > N)
    outError("Invalid parameters for hypergeometric distribution.");
  if (k > K || (n-k) > (N-K))
    return 0.0;
  double num_successes_log = binomial_coefficient_log(K, k);
  double num_failures_log = binomial_coefficient_log(N-K, n-k);
  double num_total_log = binomial_coefficient_log(N,n);
  return exp(num_successes_log + num_failures_log - num_total_log);
}

// Calculate the Frobenius norm of an N x N matrix M (flattened, rows
// concatenated) and linearly scaled by SCALE.
 double frob_norm(double m[], int n, double scale) {
   double sum = 0;
   for (int i = 0; i < n; i++) {
     for (int j = 0; j < n; j++) {
       sum += m[i*n + j] * m[i*n + j] * scale * scale;
     }
   }
   return sqrt(sum);
 }

string getOutputNameWithExt(const InputType& format, const string& output_filepath)
{
    switch (format)
    {
        case IN_MAPLE:
            return output_filepath + ".maple";
        case IN_FASTA:
            return output_filepath + ".fa";
        case IN_PHYLIP:
            return output_filepath + ".phy";
        default:
            return output_filepath + ".phy";
    }
}
