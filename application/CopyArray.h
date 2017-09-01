/*
 * CopyArray.h
 *
 *  Created on: 13-Apr-2014
 *      Author: amit
 */

#ifndef COPYARRAY_H_
#define COPYARRAY_H_

/*
 inline void CopyArrayInt(int size, int *arrSrc, int *arrDest){

 for (int i=0;i< size;i++){
 arrDest[i] = arrSrc[i];
 }

 }

 inline void CopyArrayDouble(int size, double *arrSrc, double *arrDest){

 for (int i=0;i< size;i++){
 arrDest[i] = arrSrc[i];
 }

 }
 */

/* Function ::
 * src : the Source Vector to be copied into the Vector dest vector.
 * dest : the Destination/Target Vector that contains the copy of the src vector's values.
 */
inline std::vector<double> CopyVector(std::vector<double> src) {

	std::vector<double> dest;
	int length = src.size();

	dest.resize(length);
	for (int i = 0; i < length; i++)
		dest[i] = src[i];
	return dest;
}

#endif /* COPYARRAY_H_ */
