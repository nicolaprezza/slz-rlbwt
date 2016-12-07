/*
 *  This file is part of slz-rlbwt.
 *  Copyright (c) by
 *  Nicola Prezza <nicolapr@gmail.com>,
 *  Djamal Belazzougui, Fabio Cunial, Travis Gagie, and Mathieu Raffinot
 *
 *   slz-rlbwt is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   slz-rlbwt is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 */

/*
 * slz-rlbwt.h
 *
 *  Created on: Mar 18, 2015
 *      Author: nicola
 *
 *  The slz-rlbwt data strucures combines the rlbwt and components from a LZ77 index. Space is
 *  O(r+z) words, r being the number of runs in the BWT of the text and z being the number of phrases
 *  in the LZ77 parse of the text. This version of the index also uses sparsification of LZ77 factors.
 *
 *
 */

#ifndef LZrlbwt_H_
#define LZrlbwt_H_

#include <cstdlib>
#include <iostream>
#include <vector>

#include <sdsl/wavelet_trees.hpp>
#include <rlbwt.h>
//#include <algorithms/h0_lz77.hpp>
#include <algorithms/rle_lz77_v1.hpp>
#include <sparse_sd_vector.h>
#include <range_search_2D.h>
#include <st_subset.h>
#include <packed_view.h>
#include <dynamic.hpp>
#include <stack>

using namespace std;
using namespace sdsl;
using namespace dyn;

namespace lzrlbwt{

//typedef h0_lz77<wt_fmi> lz77_t;
typedef rle_lz77_v1 lz77_t;

template<
	class fm_index_t = rlbwt_sd, 		// RLBWT using Elias-Fano compression on the gap lengths.
	class sparse_bitvector_t = sparse_sd_vector<>,		// elias-fano sparse bitvector
	class range_search_t = range_search_2D<> 	// class implementing 2-sided and 4-sided range search.
>
class slz_rlbwt{

public:

	//default empty constructor
	slz_rlbwt(){};

	// build the slz_rlbwt of the input string
	/* \param input		string to be indexed
	 * \param verbose	verbose output?
	 * \param skip:		Sparsification of LZ77: after a phrase, skip this number of characters before
	 * 					opening a new LZ77 phrase.
	 * NOTE: we assume that all letters in the input text are right-maximal.
	 */
	slz_rlbwt(string &input, bool verbose=false, ulint skip = DEFAULT_SKIP){

		this->skip = skip;

		assert(skip>0);
		assert(not contains_reserved_chars(input));

		build(input, verbose);

	}

	// return number of occurrences of input pattern
	/* \param pattern	string to be searched
	 * \return number of occurrences of the input pattern
	 */
	ulint count(string &pattern){

		range_t range = fm_index_rev.count(pattern,true);

		return range.second < range.first ? 0 : range.second - range.first + 1;

	}

	// locate all occurrences of a pattern in the indexed text. This function uses recursion, and can be space-inefficient
	// if there are a lot of occurrences
	/* \param pattern	string to be searched
	 * \return all occurrences (i.e. text positions) of the input pattern
	 *
	 * fixme: this function is bugged. It does not return correct number of occurrences
	 *
	 */
	/*vector<ulint> locate3(string &pattern, bool optimize = true){

		//first, find range of pattern (reversed because we use the reverse index)
		auto range = fm_index_rev.count(pattern,true);

		ulint n_occ = range.second < range.first ? 0 : range.second - range.first + 1;

		//we will extract text in range and put here all primary occurrences found
		vector<ulint> primary_occurrences;

		{

			vector<ulint> pocc;

			extract_occ_fwd(range, pocc, pattern.size()+skip);

			for(auto o: pocc){

				assert(o >= pattern.size() - 1);

				o = o - (pattern.size() - 1);
				primary_occurrences.push_back(o);

			}

		}

		vector<ulint> occ;

		for(auto o1:primary_occurrences){

			occ.push_back(o1);

			if(occ.size() < n_occ){

				for(auto o2 : query_2_sided_recursive(o1,o1+pattern.size()-1))
					occ.push_back(o2);

			}

		}

		assert(n_occ == occ.size());

		return occ;

	}*/

	/*vector<ulint> locate2(string &pattern, bool optimize = true){

		//first, find range of pattern (reversed because we use the reverse index)
		auto range = fm_index_rev.count(pattern,true);

		ulint n_occ = range.second < range.first ? 0 : range.second - range.first + 1;

		//we will extract text in range and put here all primary occurrences found
		vector<ulint> primary_occurrences;

		for(ulint i=range.first;i<=range.second;++i){

			auto pocc = locate_FL_until_sample_found(i, pattern.size()+skip);

			if(pocc != bwt_length){

				assert(pocc >= pattern.size() - 1);

				pocc = pocc - (pattern.size() - 1);

				primary_occurrences.push_back(pocc);

			}

		}

		vector<ulint> occ;

		for(auto o1:primary_occurrences) occ.push_back(o1);

		for(auto o1:primary_occurrences){

			if(occ.size() < n_occ){

				for(auto o2 : query_2_sided_recursive(o1,o1+pattern.size()-1))
					occ.push_back(o2);

			}

		}

		assert(n_occ == occ.size());

		return occ;

	}*/


	// locate all occurrences of a pattern in the indexed text
	// save occurrences to a stream (e.g. file) as they are found
	// internally, a stack is maintained, and the tree of occurrences copies
	// is DFS-traversed
	/* \param pattern	string to be searched
	 * \param ostr		occurrences are streamed here using 8 bytes per occurrence (uint64_t)
	 * \return all occurrences (i.e. text positions) of the input pattern
	 */
	ulint locate(string &pattern, ostream &ostr){

		auto m = pattern.length();

		//first, find range of pattern (reversed because we use the reverse index)
		auto range = fm_index_rev.count(pattern,true);

		ulint n_occ = range.second < range.first ? 0 : range.second - range.first + 1;
		ulint found_occ = 0;

		std::stack<ulint> st;//store occurrences in a stack. DFS visit of occurrences tree

		ulint p_occ = 0;	//number of primary occurrences

		for(ulint i=range.first;i<=range.second and found_occ < n_occ;++i){

			//locate if primary occurrence by extracting at most m+skip-1 characters from the end of the pattern
			auto l = locate_FL_until_sample_found(i, m+skip-1);

			if(l != bwt_length){

				assert(l >= m - 1);

				l = l - (m - 1);

			}

			//if locate succeeded
			if(l != bwt_length){

				p_occ++;

				st.push(l);

				while(not st.empty() and found_occ < n_occ){

					auto o = st.top();
					st.pop();

					//write occurrence to stream
					ostr.write((char*)&o,sizeof(o));
					found_occ++;

					{

						auto end = o+m-1;
						auto points = lz_2_sided_range_search.two_sided_range_search_ul({o,end});

						for(auto p : points){

							//p.second is a phrase rank. We need to re-map this with a select query
							assert(p.second < last_text.rank(last_text.size()));
							ulint dest_end = last_text.select(p.second);

							assert(p.first.y>=end);
							ulint shift = p.first.y-end;

							assert(dest_end >= m+shift);
							ulint j = (dest_end - m) - shift;

							//call recursively on the occurrence found
							assert(j+m-1<bwt_length-1);
							assert(j>o);
							st.push(j);

						}

					}

				}

			}

		}

		assert(n_occ == found_occ);

		return found_occ;

	}

	// save the index to file
	/* \param 	basename	basename of the output index files. Extensions will be automatically
	 * 			added to the input path
	 * \param verbose	verbose output?
	 */
	void save_to_file(string base_name, bool verbose = false){

		string fm_index_rev_basename = string(base_name).append(".lzrlbwt.fmi_rev");
		string two_sided_filename = string(base_name).append(".lzrlbwt.2_sided_range");
		string last_bwt_filename = string(base_name).append(".lzrlbwt.last_bwt");
		string last_text_filename = string(base_name).append(".lzrlbwt.last_text");
		string sa_samples_filename = string(base_name).append(".lzrlbwt.sa_samples");

		{

			if(verbose)
				cout << "Saving text positions of trailing characters of LZ factors ... " << flush;

			std::ofstream out (last_text_filename,std::ofstream::binary);
			last_text.serialize(out);
			out.close();

			if(verbose)
				cout << "done!" << endl;

		}

		{

			if(verbose)
				cout << "Saving bwt positions of trailing characters of LZ factors ... " << flush;

			std::ofstream out (last_bwt_filename,std::ofstream::binary);
			last_F.serialize(out);
			out.close();

			if(verbose)
				cout << "done!" << endl;

		}

		if(verbose)
			cout << "Saving reverse RLBWT ... " << flush;

		fm_index_rev.save_to_disk(fm_index_rev_basename);

		if(verbose)
			cout << "done!" << endl;

		if(verbose)
			cout << "Saving two-sided range data structure ... " << flush;

		lz_2_sided_range_search.save_to_file(two_sided_filename);

		if(verbose)
			cout << "done!" << endl;


		{

			if(verbose)
				cout << "Saving SA samples ... " << flush;

			std::ofstream out (sa_samples_filename,std::ofstream::binary);

			ulint sa_samples_container_size = SA_samples.container().size();
			ulint bitlength = SA_samples.width();

			out.write((char*)&sa_samples_container_size,sizeof(sa_samples_container_size));
			out.write((char*)&bitlength,sizeof(bitlength));

			out.write((char*)SA_samples.container().data(),sa_samples_container_size*sizeof(ulint));

			//save also global variables in this file
			out.write((char*)&skip,sizeof(ulint));
			out.write((char*)&terminator_pos,sizeof(ulint));
			out.write((char*)&bwt_length,sizeof(ulint));

			out.close();

			if(verbose)
				cout << "done!" << endl;

		}


	}

	// index from file
	/* \param basename_path		basename of the index files
	 * \param opt	options: verbose output / load bidirectional index (i.e. do not load reverse bwt)
	 */
	void load_from_file(string base_name, bool verbose=false){

		string fm_index_rev_basename = string(base_name).append(".lzrlbwt.fmi_rev");
		string two_sided_filename = string(base_name).append(".lzrlbwt.2_sided_range");
		string last_bwt_filename = string(base_name).append(".lzrlbwt.last_bwt");
		string last_text_filename = string(base_name).append(".lzrlbwt.last_text");
		string sa_samples_filename = string(base_name).append(".lzrlbwt.sa_samples");

		{

			if(verbose)
				cout << "Loading text positions of trailing characters of LZ factors ... " << flush;

			std::ifstream in (last_text_filename,std::ifstream::binary);
			last_text.load(in);
			in.close();

			if(verbose)
				cout << "done!" << endl;

		}

		{

			if(verbose)
				cout << "Loading bwt positions of trailing characters of LZ factors ... " << flush;

			std::ifstream in (last_bwt_filename,std::ifstream::binary);
			last_F.load(in);
			in.close();

			if(verbose)
				cout << "done!" << endl;

		}


		if(verbose)
			cout << "Loading reverse RLBWT ... " << flush;

		fm_index_rev.load_from_disk(fm_index_rev_basename);

		if(verbose)
			cout << "done!" << endl;


		if(verbose)
			cout << "Loading two-sided range data structure ... " << flush;

		lz_2_sided_range_search.load_from_file(two_sided_filename);

		if(verbose)
			cout << "done!" << endl;

		{

			if(verbose)
				cout << "Loading SA samples ... " << flush;

			std::ifstream in (sa_samples_filename,std::ifstream::binary);

			ulint sa_samples_container_size;
			ulint bitlength;

			in.read((char*)&sa_samples_container_size,sizeof(sa_samples_container_size));
			in.read((char*)&bitlength,sizeof(bitlength));

			SA_samples = packed_view<vector>(bitlength,last_F.rank(last_F.size()));

			in.read((char*)SA_samples.container().data(),sa_samples_container_size*sizeof(ulint));

			//load also global variables from this file
			in.read((char*)&skip,sizeof(ulint));
			in.read((char*)&terminator_pos,sizeof(ulint));
			in.read((char*)&bwt_length,sizeof(ulint));
			in.close();

			if(verbose)
				cout << "done!" << endl;

		}

	}

	static const ulint DEFAULT_SKIP = 64;

private:

	/*
	 * input: F range containing a unique character, a set where to put found
	 * primary occurrences, number i of FL steps left to do
	 */
	void extract_occ_fwd(range_t rn, vector<ulint>& occ, ulint i){

		if(i==0){

			for(ulint k = rn.first; k <= rn.second;++k){

				auto j = locate_fwd(k);
				if(j<bwt_length) occ.push_back(j + 1);

			}

			return;

		}

		//search and extract SA samples of marked positions inside interval rn
		find_occ(rn, occ, i);

		//map range from F to L
		vector<range_t> ranges_on_F = fm_index_rev.FL(rn);

		//recursion
		for(auto r : ranges_on_F){

			extract_occ_fwd(r,occ,i-1);

		}

	}

	/*
	 * extract skip-1 characters after this F-position
	 * returns bwt_length if no SA samples are found
	 */
	ulint locate_fwd(ulint F_pos, ulint i = 0){

		if(i==skip-1) return bwt_length;

		if(last_F.at(F_pos)){

			return last_text.select(SA_samples[last_F.rank(F_pos)]) + i;

		}

		return locate_fwd(fm_index_rev.FL(F_pos),i+1);

	}

	/*
	 * extract i characters after this F-position. Stop when a sample is found.
	 * returns bwt_length if no SA samples are found
	 */
	ulint locate_FL_until_sample_found(ulint F_pos, ulint i = 0){

		if(i==0) return bwt_length;

		if(last_F.at(F_pos)){//sample found

			return last_text.select(SA_samples[last_F.rank(F_pos)]);

		}


		//sample not found: search forward in the rev BWT (i.e. backward in the text)
		auto l = locate_FL_until_sample_found(fm_index_rev.FL(F_pos),i-1);

		return l == bwt_length ? bwt_length : l + 1;

	}

	/*
	 * extract i characters before this F-position. Stop when a sample is found.
	 * returns bwt_length if no SA samples are found
	 */
	ulint locate_LF_until_sample_found(ulint F_pos, ulint i = 0){

		if(i==0) return bwt_length;

		if(last_F.at(F_pos)){//sample found

			return last_text.select(SA_samples[last_F.rank(F_pos)]);

		}


		//sample not found: search backward
		auto l = locate_LF_until_sample_found(fm_index_rev.LF(F_pos),i-1);

		return l == bwt_length ? bwt_length : l - 1;

	}

	void find_occ(range_t rn, vector<ulint>& occ, ulint i){

		assert(i>0);
		assert(rn.first>0); //# must not be contained in the range

		assert(rn.first<=last_F.size());
		assert(rn.second+1<=last_F.size());
		assert(SA_samples.size() == last_F.rank(last_F.size()));
		assert(SA_samples.size() == last_text.rank(last_text.size()));

		for(ulint j = last_F.rank(rn.first);j<last_F.rank(rn.second+1);++j){

			assert(j<SA_samples.size());
			ulint text_position = last_text.select(SA_samples[j]);

			assert(text_position >= (i - 1));
			occ.push_back( text_position - (i - 1) );

		}

	}

	void build(string &input, bool verbose){

		//Append LZ77 terminator
		//0 and 1 are reserved by the RLBWT and sdsl, so
		//we use 2 for the LZ terminator
		input = input + char(2);

		//size of the text
		ulint text_length = input.size();
		bwt_length=text_length + 1;

		if(verbose)
			cout << "Building reverse RLBWT ... " << flush;

		{
			string rev_text;

			for(ulint i=0;i<text_length;++i)
				rev_text.push_back(input[text_length - i - 1]);

			fm_index_rev = fm_index_t(rev_text);
		}

		if(verbose)
			cout << "done!" << endl;

		if(verbose)
			cout << "Initializing structures for the LZ77 parser ... " << endl;

		lz77_t parser;

		{

			std::stringstream iss(input);

			//detect characters probabilities (for Huffman compression)
			//parser = lz77_t(iss,32);

			parser = lz77_t(iss);

			string rev_bwt = fm_index_rev.toString();

			parser.load_bwt(rev_bwt, fm_index_rev.get_terminator(), verbose);

		}



		{

			//std::stringstream iss(input);
			std::ostringstream oss;

			parser.bwt_to_lz77(oss,skip,verbose);

			//parser.parse(iss,oss,skip,true);

			//this vector will contain true on trailing characters (on text)
			vector<bool> last_text_vec;

			ulint z=0;	//number of LZ phrases

			//this vector will contain the 2D points for range search, one per phrase
			vector<pair<point_2d_t, ulint> > two_sided_points;

			assert(oss.str().length() % 17 == 0);//8+8+1 bytes per factor

			//auto phrases = oss.str().c_str();
			ulint n_phrases = oss.str().length()/17;

			istringstream iss2(oss.str());
			//oss.seekp(0, std::ios::beg);

			//ulint i = 0;

			//scan parse
			while(z<n_phrases){

				ulint start_pos;
				ulint len;
				char c;

				iss2.read((char *)&start_pos, 8);
				iss2.read((char *)&len, 8);
				iss2.read(&c, 1);

				//cout << " " << start_pos << " " << len << " " << c << endl;

				assert(start_pos < text_length);
				assert(len < text_length);

				//append false for each copied character in the phrase
				for(ulint j=0;j<len;++j)
					last_text_vec.push_back(false);

				//append a true: this is the first skipped character after the phrase
				last_text_vec.push_back(true);

				//append false for each skipped character. Stop if we go beyond
				//text length
				for(ulint j=1;j<skip and last_text_vec.size() < text_length;++j)
					last_text_vec.push_back(false);

				//for all factors with a start position, insert the corresponding point in the
				//2-sided range search structure
				if(len>0){

					//pair: <start_position, end_position> of the copied string
					point_2d_t point = {start_pos, start_pos + len - 1};

					//insert the point in the vector of points.
					//we associate phrase rank (z) to each point
					two_sided_points.push_back({point,z});

				}

				z++;	//increment number of phrases

			}

			if(verbose)
				cout << "Building two-sided range data structure ... " << flush;

			//build range search structure. Do not re-map y coordinates: saves z words of space
			lz_2_sided_range_search = range_search_t(two_sided_points,true,false);

			if(verbose)
				cout << "done!" << endl;

			assert(last_text_vec.size()==text_length);

			last_text = sparse_bitvector_t(last_text_vec);

		}

		// build array 'last'

		{

			if(verbose)
				cout << "Building sparse bitvector marking end of phrases on BWT ... " << flush;

			//first, build a standard bitvector where phrase ends (on the F column of the BWT of the
			//reversed text) are marked with a 1. Then, convert it to a sparse bitvector.
			vector<bool> last_F_vec(bwt_length, false);

			//position on column F
			ulint F_pos = fm_index_rev.LF(0);

			//navigate the reverse BWT, reading the (forward) text from first to last_F character
			for(ulint i=0;i<text_length;++i){

				if(last_text[i]){

					assert(F_pos < last_F_vec.size());
					last_F_vec[F_pos] = true;		//mark the position

				}

				F_pos = fm_index_rev.LF(F_pos);		//update position in the BWT

			}

			//now bwt_pos is the position on F of #
			assert(F_pos==0);

			//now convert the bitvector to a sparse bitvector
			last_F = sparse_bitvector_t(last_F_vec);

			if(verbose)
				cout << "done!" << endl;

		}

		{

			if(verbose)
				cout << "Sampling suffix array at the end of LZ phrases ... " << flush;

			assert(text_length>0);

			//in SA_samples we write phrase ranks (on the text)
			ulint z = last_text.rank(last_text.size());

			assert( z == last_F.rank(last_F.size()));

			//number of bits required to write the number z
			ulint bitlength =  64 - __builtin_clzll(z);

			SA_samples = packed_view<vector>(bitlength, z);
			assert(SA_samples.width()==bitlength);

			//position on F column of rev BWT
			ulint F_pos = fm_index_rev.LF(0);

			//navigate the reverse BWT, reading the (forward) text from first to last character
			//i = text position corresponding to F_pos
			for(ulint i=0;i<text_length;++i){

				//sampled position
				assert(F_pos < last_F.size());
				if(last_F.at(F_pos)){

					//to save bits, we actually sample the rank of current phrase (+1 because we are
					//after the bit set of this phrase beginning in 'begin_of_phrase')
					SA_samples[last_F.rank(F_pos)] = last_text.rank(i);

				}

				F_pos = fm_index_rev.LF(F_pos);		//update position in the BWT

			}

			//now bwt_pos is the position on F of #
			assert(F_pos==0);

			if(verbose)
				cout << "done!" << endl;

		}

	}


	/*
	 * find secondary occurrences copied from T[begin,...,end].
	 * Calls recursively itself until no more occurrences can be found
	 */
	vector<ulint> query_2_sided_recursive(ulint begin, ulint end){

		vector<ulint> result;

		assert(end>=begin);

		//pattern length
		ulint m = end-begin+1;

		//launch a 2 sided range search
		//result: a vector of pairs < <b,e>, i >, meaning that T[b,e] is copied in T[i,...,i+(e-b)]

		assert(end<bwt_length-1);

		auto points = lz_2_sided_range_search.two_sided_range_search_ul({begin,end});

		for(auto p : points){

			//p.second is a phrase rank. We need to re-map this with a select query
			assert(p.second < last_text.rank(last_text.size()));
			ulint dest_end = last_text.select(p.second);

			assert(p.first.y>=end);
			ulint shift = p.first.y-end;

			assert(dest_end >= m+shift);
			ulint j = (dest_end - m) - shift;

			//call recursively on the occurrence found
			assert(j+m-1<bwt_length-1);
			assert(j>begin);
			auto recursive_occ = query_2_sided_recursive(j,j+m-1);

			//append all found occurrences to the final result
			result.push_back(j);

			for(auto x : recursive_occ)
				result.push_back(x);

		}

		return result;

	}

	bool contains_reserved_chars(string &s){

		for(auto c:s) if(c==0 or c==1 or c==2) return true;

		return false;

	}

	ulint bwt_length = 0;		//length of the BWT (i.e. 1+text length)
	ulint terminator_pos = 0;	//position of terminator character in the forward BWT

	range_search_t lz_2_sided_range_search; 	// points with occurring positions of factors

	fm_index_t fm_index_rev;	// O(R) words. FM index for the reversed text

	//z words. marks with a 1 positions in the F column of rev BWT
	//that correspond to the first character after a phrase
	sparse_bitvector_t last_F;
	//z words. marks with a 1 positions in the TEXT
	//that correspond to the first character after a phrase. Each bit
	//set here corresponds to a bit set in the array last_F
	sparse_bitvector_t last_text;

	//used only if index is of light type (otherwise, samples
	//are stored inside lz_4_sided_range_search).
	//stores ranks of bits set in last_text
	packed_view<vector> SA_samples;

	ulint skip = DEFAULT_SKIP;	//sparsification of LZ77: after a phrase, skip this number
								//of characters (including trailing character), then begin new phrase

};//class slz_rlbwt

}//namespace lzrlbwt

#endif /* LZrlbwt_H_ */
