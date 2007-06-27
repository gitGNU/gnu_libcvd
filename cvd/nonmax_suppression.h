#ifndef CVD_NONMAX_SUPPRESSION_H
#define CVD_NONMAX_SUPPRESSION_H

#include <vector>
#include <cvd/image_ref.h>

namespace CVD
{
	/**Perform nonmaximal suppression on a set of features.

	@param corners The corner locations
	@param scores  The corners' scores
	@param max_corners The locally maximal corners.
	@ingroup gVision
	*/
	void nonmax_suppression(const vector<ImageRef>& corners, const vector<int>& scores, vector<ImageRef>& nmax_corners);


	/**Perform nonmaximal suppression on a set of features.

	@param corners The corner locations
	@param scores  The corners' scores
	@param max_corners The locally maximal corners, and their scores.
	@ingroup gVision
	*/
	void nonmax_suppression_with_scores(const vector<ImageRef>& corners, const vector<int>& socres, vector<pair<ImageRef,int> >& max_corners);

}
