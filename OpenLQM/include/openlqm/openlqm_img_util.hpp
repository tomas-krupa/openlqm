/*
* Copyright 2025 Noblis, Inc.
* https://fingerprint.nist.gov/openlqm
* 
* Licensed under the Apache License, Version 2.0 (the "License"); you may not
* use this file except in compliance with the License. You may obtain a copy of
* the License at https://www.apache.org/licenses/LICENSE-2.0.
* 
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
* License for the specific language governing permissions and limitations under
* the License.
* 
* OpenLQM was developed by Noblis, Inc. under contract to the National
* Institute of Standards and Technology in 2024-2025, as a port of "LQMetric,"
* which was developed by Noblis, Inc. under contract to the Federal Bureau of
* Investigation's Criminal Justice Information Services Division in 2012-2014.
*/

#pragma once

#include "openlqm.hpp"
#include "openlqm_minutia.hpp"
#include <vector>
#include <utility>

// Initial value for uncomputed Direction Map cells
#define LQM_INVALID_DIR -1

// Integer values corresponding to `true` and `false` in various contexts
#define LQM_TRUE 1
#define LQM_FALSE 0

// Integer values returned by some functions that scan pixel neighbors
#define LQM_FOUND                 LQM_TRUE
#define LQM_NOT_FOUND             LQM_FALSE

// Return codes representing the relative positioning between a chosen pair of pixels
#define LQM_NORTH                0
#define LQM_SOUTH                4
#define LQM_EAST                 2
#define LQM_WEST                 6

// Maximum value of an unsigned 6-bit gray pixel
#define LQM_MAX_6_BIT_GRAY 64

/* Minimum total DFT power for any given block  */
/* which is used to compute an average power.   */
/* By setting a non-zero minimum total,possible */
/* division by zero is avoided.                 */
#define LQM_MIN_POWER_SUM           10.0

namespace OpenLQM {
   namespace Core {
      /* Determine the maximum amount of image padding required to support */
      /* LQM processes.                                                    */
      int GetMaxPadding(const int map_windowsize, const int map_windowoffset, const int dirbin_grid_w, const int dirbin_grid_h);

      /*************************************************************************
      **************************************************************************
      PadUCharImage - Copies an 8-bit grayscale images into a larger
                        output image centering the input image so as to
                        add a specified amount of pixel padding along the
                        entire perimeter of the input image.  The amount of
                        pixel padding and the intensity of the pixel padding
                        are specified.  An alternative to padding with a
                        constant intensity would be to copy the edge pixels
                        of the centered image into the adjacent pad area.

         Input:
            idata     - input 8-bit grayscale image
            iw        - width (in pixels) of the input image
            ih        - height (in pixels) of the input image
            pad       - size of padding (in pixels) to be added
            pad_value - intensity of the padded area
         Output:
            out_paddedData    - points to the newly padded image
            ow                - width (in pixels) of the padded image
            oh                - height (in pixels) of the padded image
      **************************************************************************/
      void PadUCharImage(std::vector<unsigned char>& out_paddedData, int& ow, int& oh,
                        unsigned char *idata, const int iw, const int ih,
                        const int pad, const int pad_value);

      /*************************************************************************
      **************************************************************************
      Bits8to6 - Takes an array of unsigned characters and bitwise shifts
                  each value 2 positions to the right.  This is equivalent
                  to dividing each value by 4.  This puts original values
                  on the range [0..256) now on the range [0..64).  Another
                  way to say this, is the original 8-bit values now fit in
                  6 bits.  I would really like to make this dependency
                  go away.

         Input:
            idata - input array of unsigned characters
            iw    - width (in characters) of the input array
            ih    - height (in characters) of the input array
         Output:
            idata - contains the bit-shifted results
      **************************************************************************/
      void Bits8to6(unsigned char *idata, const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      GenImageMaps - Computes a set of image maps. The first map is a
               Direction Map which is a 2D vector of integer directions,
               where each direction represents the dominant ridge flow
               in a block of the input grayscale image. The Low Contrast
               Map flags blocks with insufficient contrast. The Low
               Flow Map flags blocks with insufficient ridge flow.
               The High Curve Map flags blocks containing high curvature.
               This routine will generate maps for an arbitrarily sized,
               non-square, image.

         Input:
            pdata     - padded input image data (8 bits [0..256) grayscale)
            pw        - padded width (in pixels) of the input image
            ph        - padded height (in pixels) of the input image
            dir2rad   - lookup table for converting integer directions
            dftwaves  - structure containing the DFT wave forms
            dftgrids  - structure containing the rotated pixel grid offsets
            lqmParams - parameters and thresholds for controlling OpenLQM
         Output:
            out_dmap     - points to the created Direction Map
            out_lcmap    - points to the created Low Contrast Map
            out_lfmap    - points to the Low Ridge Flow Map
            out_hcmap    - points to the High Curvature Map
            omw       - width (in blocks) of the maps
            omh       - height (in blocks) of the maps
      **************************************************************************/
      void GenImageMaps(std::vector<int>& out_dmap, std::vector<int>& out_lcmap, std::vector<int>& out_lfmap, std::vector<int>& out_hcmap,
                  int *omw, int *omh,
                  unsigned char *pdata, const int pw, const int ph,
                  const Dir2Rad& dir2rad, const DFTWaves& dftwaves,
                  const RotGrids& dftgrids, const LQMParams& lqmParams,
                  unsigned char *gray10pctMap, unsigned char *gray90pctMap, unsigned char *grayMedianMap, unsigned char *grayRangeMap, unsigned char *grayCountMap,
                  double *maxMagnitude, double *normMagnitude, double *lowFreqMagnitude,
                  unsigned char *directionChangeMap,
                  unsigned char *validNeighbors, unsigned char *curvatureMap);

      /*************************************************************************
      **************************************************************************
      BlockOffsets - Divides an image into mw X mh equally sized blocks,
            returning a list of offsets to the top left corner of each block.
            For images that are even multiples of BLOCKSIZE, blocks do not
            not overlap and are immediately adjacent to each other.  For image
            that are NOT even multiples of BLOCKSIZE, blocks continue to be
            non-overlapping up to the last column and/or last row of blocks.
            In these cases the blocks are adjacent to the edge of the image and
            extend inwards BLOCKSIZE units, overlapping the neighboring column
            or row of blocks.  This routine also accounts for image padding
            which makes things a little more "messy". This routine is primarily
            responsible for providing the ability to processs arbitrarily-sized
            images.  The strategy used here is simple, but others are possible.

         Input:
            iw                - width (in pixels) of the orginal input image
            ih                - height (in pixels) of the orginal input image
            pad               - the padding (in pixels) required to support the desired
                              range of block orientations for DFT analysis.  This padding
                              is required along the entire perimeter of the input image.
                              For certain applications, the pad may be zero.
            blocksize         - the width and height (in pixels) of each image block
         Output:
            out_blockOffsets  - the list of pixel offsets to the origin of
                              each block in the "padded" input image
            ow                - the number of horizontal blocks in the input image
            oh                - the number of vertical blocks in the input image
         Return Code:
            Zero              - successful completion
            Negative          - system error
      **************************************************************************/
      void BlockOffsets(std::vector<int>& out_blockOffsets, int& ow, int& oh,
               const int iw, const int ih, const int pad, const int blocksize);


      /*************************************************************************
      **************************************************************************
      GenInitialMaps - Creates an initial Direction Map from the given
                  input image.  It very important that the image be properly
                  padded so that rotated grids along the boundary of the image
                  do not access unkown memory.  The rotated grids are used by a
                  DFT-based analysis to determine the integer directions
                  in the map. Typically this initial vector of directions will
                  subsequently have weak or inconsistent directions removed
                  followed by a smoothing process.  The resulting Direction
                  Map contains valid directions >= 0 and INVALID values = -1.
                  This routine also computes and returns 2 other image maps.
                  The Low Contrast Map flags blocks in the image with
                  insufficient contrast.  Blocks with low contrast have a
                  corresponding direction of INVALID in the Direction Map.
                  The Low Flow Map flags blocks in which the DFT analyses
                  could not determine a significant ridge flow.  Blocks with
                  low ridge flow also have a corresponding direction of
                  INVALID in the Direction Map.

         Input:
            blkoffs   - offsets to the pixel origin of each block in the padded image
            mw        - number of blocks horizontally in the padded input image
            mh        - number of blocks vertically in the padded input image
            pdata     - padded input image data (8 bits [0..256) grayscale)
            pw        - width (in pixels) of the padded input image
            ph        - height (in pixels) of the padded input image
            dftwaves  - structure containing the DFT wave forms
            dftgrids  - structure containing the rotated pixel grid offsets
            lqmParams - parameters and thresholds for controlling OpenLQM
         Output:
            out_dmap     - points to the newly created Direction Map
            out_lcmap    - points to the newly created Low Contrast Map
            out_lfmap    - points to the newly created Low Flow Map
         Return Code:
            Zero     - successful completion
            Negative - system error
      **************************************************************************/
      void GenInitialMaps(std::vector<int>& out_dmap, std::vector<int>& out_lcmap, std::vector<int>& out_lfmap,
                     int *blkoffs, const int mw, const int mh,
                     unsigned char *pdata, const int pw, const int ph,
                     const DFTWaves& dftWaves, const RotGrids& dftGrids,
                     const LQMParams& lqmParams,
                     unsigned char *gray10pctMap, unsigned char *gray90pctMap, unsigned char *grayMedianMap, unsigned char *grayRangeMap, unsigned char *grayCountMap,
                     double *maxMagnitude, double *normMagnitude, double *lowFreqMagnitude);

      /*************************************************************************
      **************************************************************************
      AllocDirPowers - Allocates the memory associated with DFT power
               vectors.  The DFT analysis is conducted block by block in the
               input image, and within each block, N wave forms are applied
               at M different directions.

         Input:
            nwaves - number of DFT wave forms
            ndirs  - number of orientations (directions) used in DFT analysis
         Output:
            out_powers - pointer to the allocated power vectors
      **************************************************************************/
      void AllocDirPowers(std::vector<std::vector<double>>& out_powers, const int nwaves, const int ndirs);

      /*************************************************************************
      **************************************************************************
      AllocPowerStats - Allocates memory associated with set of statistics
                  derived from DFT power vectors computed in a block of the
                  input image.  Statistics are not computed for the lowest DFT
                  wave form, so the length of the statistics arrays is 1 less
                  than the number of DFT wave forms used.  The staistics
                  include the Maximum power for each wave form, the direction
                  at which the maximum power occured, and a normalized value
                  for the maximum power.  In addition, the statistics are
                  ranked in descending order based on normalized squared
                  maximum power.

         Input:
            nstats - the number of waves forms from which statistics are to be
                     derived (N Waves - 1)
         Output:
            out_wis           - output array holding the ranked wave form indices
                              of the corresponding statistics
            out_powmaxs       - output array holding the maximum DFT power for each
                              wave form
            out_powmaxDirs    - output array holding the direction corresponding to
                              each maximum power value
            out_pownorms      - output array holding the normalized maximum power
      **************************************************************************/
      void AllocPowerStats(std::vector<int>& out_wis, std::vector<double>& out_powmaxs, std::vector<int>& out_powmaxDirs,
                           std::vector<double>& out_pownorms, const int nstats);

      /*************************************************************************
      LowContrastBlock - Takes the offset to an image block of specified
                  dimension, and analyzes the pixel intensities in the block
                  to determine if there is sufficient contrast for further
                  processing.

         Input:
            blkoffset - byte offset into the padded input image to the origin of
                        the block to be analyzed
            blocksize - dimension (in pixels) of the width and height of the block
                        (passing separate blocksize from LQMParams on purpose)
            pdata     - padded input image data (8 bits [0..256) grayscale)
            pw        - width (in pixels) of the padded input image
            lqmParams  - parameters and thresholds for controlling LQM
         Output:
            pct10      - 10th percentile gray level (6-bit)
            pct90      - 90th percentile gray level (6-bit)
            median     - median gray level (6-bit)
            range      - range of gray levels (6-bit)
            numgrays   - count of distinct gray values (6-bit)
         Return Value:
            true     - block has sufficiently low contrast
            false    - block has sufficiently hight contrast
      **************************************************************************
      **************************************************************************/
      bool LowContrastBlock(const int blkoffset, const int blocksize,
                           unsigned char *pdata, const int pw,
                           const LQMParams& lqmParams,
                           unsigned char *pct10, unsigned char *pct90, unsigned char *median, unsigned char *range, unsigned char *numgrays);

      /*************************************************************************
      **************************************************************************
      DFTDirPowers - Conducts the DFT analysis on a block of image data.
            The image block is sampled across a range of orientations
            (directions) and multiple wave forms of varying frequency are
            applied at each orientation.  At each orentation, pixels are
            accumulated along each rotated pixel row, creating a vector
            of pixel row sums.  Each DFT wave form is then applied
            individually to this vector of pixel row sums.  A DFT power
            value is computed for each wave form (frequency0 at each
            orientaion within the image block.  Therefore, the resulting DFT
            power vectors are of dimension (N Waves X M Directions).
            The power signatures derived form this process are used to
            determine dominant direction flow within the image block.

         Input:
            pdata     - the padded input image.  It is important that the image
                        be properly padded, or else the sampling at various block
                        orientations may result in accessing unkown memory.
            blkoffset - the pixel offset form the origin of the padded image to
                        the origin of the current block in the image
            dftWaves  - structure containing the DFT wave forms
            dftGrids  - structure containing the rotated pixel grid offsets
         Output:
            powers    - DFT power computed from each wave form frequencies at each
                        orientation (direction) in the current image block
      **************************************************************************/
      void DFTDirPowers(std::vector<std::vector<double>>& powers, unsigned char *pdata,
                     const int blkoffset,
                     const DFTWaves& dftWaves, const RotGrids& dftGrids);

      /*************************************************************************
      **************************************************************************
      SumRotBlockRows - Computes a vector or pixel row sums by sampling
                  the current image block at a given orientation.  The
                  sampling is conducted using a precomputed set of rotated
                  pixel offsets (called a grid) relative to the orgin of
                  the image block.

         Input:
            blkptr    - the pixel address of the origin of the current image block
            grid_offsets - the rotated pixel offsets for a block-sized grid
                        rotated according to a specific orientation
            blocksize - the width and height of the image block and thus the size
                        of the rotated grid
         Output:
            rowsums   - the resulting vector of pixel row sums
      **************************************************************************/
      void SumRotBlockRows(std::vector<int>& rowsums, const unsigned char *blkptr,
                              const int *grid_offsets, const int blocksize);

      /*************************************************************************
      **************************************************************************
      DFTPower - Computes the DFT power by applying a specific wave form
                  frequency to a vector of pixel row sums computed from a
                  specific orientation of the block image

         Input:
            rowsums - accumulated rows of pixels from within a rotated grid
                     overlaying an input image block
            wave    - the wave form (cosine and sine components) at a specific
                     frequency
            wavelen - the length of the wave form (must match the height of the
                     image block which is the length of the rowsum vector)
         Output:
            power   - the computed DFT power for the given wave form at the
                     given orientation within the image block
      **************************************************************************/
      void DFTPower(double *power, const int *rowsums,
                     const DFTWave& wave, const int wavelen);

      /*************************************************************************
      **************************************************************************
      DFTPowerStats - Derives statistics from a set of DFT power vectors.
               Statistics are computed for all but the lowest frequency
               wave form, including the Maximum power for each wave form,
               the direction at which the maximum power occured, and a
               normalized value for the maximum power.  In addition, the
               statistics are ranked in descending order based on normalized
               squared maximum power.  These statistics are fundamental
               to selecting a dominant direction flow for the current
               input image block.

         Input:
            powers   - DFT power vectors (N Waves X M Directions) computed for
                     the current image block from which the values in the
                     statistics arrays are derived
            fw       - the beginning of the range of wave form indices from which
                     the statistcs are to derived
            tw       - the ending of the range of wave form indices from which
                     the statistcs are to derived (last index is tw-1)
            ndirs    - number of orientations (directions) at which the DFT
                     analysis was conducted
         Output:
            wis      - list of ranked wave form indicies of the corresponding
                     statistics based on normalized squared maximum power. These
                     indices will be used as indirect addresses when processing
                     the power statistics in descending order of "dominance"
            powmaxs  - array holding the maximum DFT power for each wave form
                     (other than the lowest frequecy)
            powmax_dirs - array to holding the direction corresponding to
                        each maximum power value in powmaxs
            pownorms - array to holding the normalized maximum powers corresponding
                     to each value in powmaxs
      **************************************************************************/
      void DFTPowerStats(int *wis, double *powmaxs, int *powmax_dirs,
                           double *pownorms, std::vector<std::vector<double>>& powers,
                           const int fw, const int tw, const int ndirs);

      /*************************************************************************
      **************************************************************************
      GetMaxNorm - Analyses a DFT power vector for a specific wave form
                     applied at different orientations (directions) to the
                     current image block.  The routine retuns the maximum
                     power value in the vector, the direction at which the
                     maximum occurs, and a normalized power value.  The
                     normalized power is computed as the maximum power divided
                     by the average power across all the directions.  These
                     simple statistics are fundamental to the selection of
                     a dominant direction flow for the image block.

         Input:
            power_vector - the DFT power values derived form a specific wave form
                           applied at different directions
            ndirs      - the number of directions to which the wave form was applied
         Output:
            powmax     - the maximum power value in the DFT power vector
            powmax_dir - the direciton at which the maximum power value occured
            pownorm    - the normalized power corresponding to the maximum power
      **************************************************************************/
      void GetMaxNorm(double *powmax, int *powmax_dir,
                     double *pownorm, const double *power_vector, const int ndirs);

      /*************************************************************************
      **************************************************************************
      SortDFTWaves - Creates a ranked list of DFT wave form statistics
                     by sorting on the normalized squared maximum power.

         Input:
            powmaxs  - maximum DFT power for each wave form used to derive
                     statistics
            pownorms - normalized maximum power corresponding to values in powmaxs
            nstats   - number of wave forms used to derive statistics (N Wave - 1)
         Output:
            wis      - sorted list of indices corresponding to the ranked set of
                     wave form statistics.  These indices will be used as
                     indirect addresses when processing the power statistics
                     in descending order of "dominance"
      **************************************************************************/
      void SortDFTWaves(int *wis, const double *powmaxs, const double *pownorms,
                        const int nstats);

      /***************************************************************************
      **************************************************************************
      BubbleSortDoubleDec - Conducts a simple bubble sort returning a list
            of ranks in decreasing order and their associated items in sorted
            order as well.

         Input:
            ranks - list of values to be sorted
            items - list of items, each corresponding to a particular rank value
            len   - length of the lists to be sorted
         Output:
            ranks - list of values sorted in descending order
            items - list of items in the corresponding sorted order of the ranks.
                  If these items are indices, upon return, they may be used as
                  indirect addresses reflecting the sorted order of the ranks.
      ****************************************************************************/
      void BubbleSortDoubleDec(double *ranks, int *items,  const int len);

      /*************************************************************************
      **************************************************************************
      PrimaryDirTest - Applies the primary set of criteria for selecting
                        an IMAP integer direction from a set of DFT results
                        computed from a block of image data

         Input:
            powers      - DFT power computed from each (N) wave frequencies at each
                        rotation direction in the current image block
            wis         - sorted order of the highest N-1 frequency power statistics
            powmaxs     - maximum power for each of the highest N-1 frequencies
            powmax_dirs - directions associated with each of the N-1 maximum powers
            pownorms    - normalized power for each of the highest N-1 frequencies
            nstats      - N-1 wave frequencies (where N is the length of dft_coefs)
            lqmParams    - parameters and thresholds for controlling LQM
         Return Code:
            Zero or Positive - The selected IMAP integer direction
            LQM_INVALID_DIR - IMAP Integer direction could not be determined
      **************************************************************************/
      int PrimaryDirTest(std::vector<std::vector<double>>& powers, const int *wis,
                  const double *powmaxs, const int *powmax_dirs,
                  const double *pownorms, const int nstats,
                  const LQMParams& lqmParams,
                  double *maxMagnitude, double *normMagnitude, double *lowFreqMagnitude);

      /*************************************************************************
      **************************************************************************
      SecondaryForkTest - Applies a secondary set of criteria for selecting
                        an IMAP integer direction from a set of DFT results
                        computed from a block of image data.  This test
                        analyzes the strongest power statistics associated
                        with a given frequency and direction and analyses
                        small changes in direction to the left and right to
                        determine if the block contains a "fork".

         Input:
            powers      - DFT power computed from each (N) wave frequencies at each
                        rotation direction in the current image block
            wis         - sorted order of the highest N-1 frequency power statistics
            powmaxs     - maximum power for each of the highest N-1 frequencies
            powmax_dirs - directions associated with each of the N-1 maximum powers
            pownorms    - normalized power for each of the highest N-1 frequencies
            lqmParams    - parameters and thresholds for controlling LQM
         Return Code:
            Zero or Positive - The selected IMAP integer direction
            INVALID_DIR - IMAP Integer direction could not be determined
      **************************************************************************/
      int SecondaryForkTest(std::vector<std::vector<double>>& powers, const int *wis,
                  const double *powmaxs, const int *powmax_dirs,
                  const double *pownorms,
                  const LQMParams& lqmParams,
                  double *maxMagnitude, double *normMagnitude, double *lowFreqMagnitude);

      /*************************************************************************
      **************************************************************************
      MorphTFMap - Takes a 2D vector of TRUE and FALSE values integers
                  and dialates and erodes the map in an attempt to fill
                  in voids in the map.

         Input:
            tfmap    - vector of integer block values
            mw       - width (in blocks) of the map
            mh       - height (in blocks) of the map
         Output:
            tfmap    - resulting morphed map
      **************************************************************************/
      void MorphTFMap(int *tfmap, const int mw, const int mh);

      /*************************************************************************
      **************************************************************************
      DilateUCharImage - Dilates an 8-bit image by setting false pixels to
                  one if any of their 4 neighbors is non-zero.  Allocation
                  of the output image is the responsibility of the caller.
                  The input image remains unchanged.

         Input:
            inp       - input 8-bit image to be dilated
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
         Output:
            out       - contains to the resulting dilated image
      **************************************************************************/
      void DilateUCharImage(unsigned char *inp, unsigned char *out,
                           const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      ErodeUCharImage - Erodes an 8-bit image by setting true pixels to zero
                  if any of their 4 neighbors is zero.  Allocation of the
                  output image is the responsibility of the caller.  The
                  input image remains unchanged.  This routine will NOT
                  erode pixels indiscriminately along the image border.

         Input:
            inp       - input 8-bit image to be eroded
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
         Output:
            out       - contains the resulting eroded image
      **************************************************************************/
      void ErodeUCharImage(unsigned char *inp, unsigned char *out,
                           const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      GetSouthUChar - Returns the value of the 8-bit image pixel 1 below the
                     current pixel if defined else it returns defaultVal.

         Input:
            ptr       - points to current pixel in image
            row       - y-coord of current pixel
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            defaultVal - value returned for out-of-bounds pixel
         Return Code:
            Zero      - if neighboring pixel is undefined
                        (outside of image boundaries)
            Pixel     - otherwise, value of neighboring pixel
      **************************************************************************/
      unsigned char GetSouthUChar(unsigned char *ptr, const int row, const int iw, const int ih, unsigned char defaultVal);

      /*************************************************************************
      **************************************************************************
      GetNorthUChar - Returns the value of the 8-bit image pixel 1 above the
                     current pixel if defined else it returns defaultVal.

         Input:
            ptr       - points to current pixel in image
            row       - y-coord of current pixel
            iw        - width (in pixels) of image
            defaultVal - value returned for out-of-bounds pixel
         Return Code:
            Zero      - if neighboring pixel is undefined
                        (outside of image boundaries)
            Pixel     - otherwise, value of neighboring pixel
      **************************************************************************/
      unsigned char GetNorthUChar(unsigned char *ptr, const int row, const int iw, unsigned char defaultVal);

      /*************************************************************************
      **************************************************************************
      GetEastUChar - Returns the value of the 8-bit image pixel 1 right of the
                  current pixel if defined else it returns defaultVal.

         Input:
            ptr       - points to current pixel in image
            col       - x-coord of current pixel
            iw        - width (in pixels) of image
            defaultVal - value returned for out-of-bounds pixel
         Return Code:
            Zero      - if neighboring pixel is undefined
                        (outside of image boundaries)
            Pixel     - otherwise, value of neighboring pixel
      **************************************************************************/
      unsigned char GetEastUChar(unsigned char *ptr, const int col, const int iw, unsigned char defaultVal);

      /*************************************************************************
      **************************************************************************
      GetWestUChar - Returns the value of the 8-bit image pixel 1 left of the
                  current pixel if defined else it returns defaultVal.

         Input:
            ptr       - points to current pixel in image
            col       - x-coord of current pixel
            defaultVal - value returned for out-of-bounds pixel
         Return Code:
            Zero      - if neighboring pixel is undefined
                        (outside of image boundaries)
            Pixel     - otherwise, value of neighboring pixel
      **************************************************************************/
      unsigned char GetWestUChar(unsigned char *ptr, const int col, unsigned char defaultVal);

      /*************************************************************************
      **************************************************************************
      RemoveInconDirs - Takes a vector of integer directions and removes
                  individual directions that are too weak or inconsistent.
                  Directions are tested from the center of the IMAP working
                  outward in concentric squares, and the process resets to
                  the center and continues until no changes take place during
                  a complete pass.

         Input:
            imap      - vector of IMAP integer directions
            mw        - width (in blocks) of the IMAP
            mh        - height (in blocks) of the IMAP
            dir2rad   - lookup table for converting integer directions
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            imap      - vector of pruned input values
      **************************************************************************/
      void RemoveInconDirs(int *imap, const int mw, const int mh,
                  const Dir2Rad& dir2rad, const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      RemoveDir - Determines if an IMAP direction should be removed based
                  on analyzing its adjacent neighbors

         Input:
            imap      - vector of IMAP integer directions
            mx        - IMAP X-coord of the current direction being tested
            my        - IMPA Y-coord of the current direction being tested
            mw        - width (in blocks) of the IMAP
            mh        - height (in blocks) of the IMAP
            dir2rad   - lookup table for converting integer directions
            lqmParams - parameters and thresholds for controlling LQM
         Return Code:
            Positive - direction should be removed from IMAP
            Zero     - direction should NOT be remove from IMAP
      **************************************************************************/
      int RemoveDir(int *imap, const int mx, const int my,
                     const int mw, const int mh, const Dir2Rad& dir2rad,
                     const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      Average8NeighDir - Given an IMAP direction, computes an average
                        direction from its adjacent 8 neighbors returning
                        the average direction, its strength, and the
                        number of valid direction in the neighborhood.

         Input:
            imap      - vector of IMAP integer directions
            mx        - IMAP X-coord of the current direction
            my        - IMPA Y-coord of the current direction
            mw        - width (in blocks) of the IMAP
            mh        - height (in blocks) of the IMAP
            dir2rad   - lookup table for converting integer directions
         Output:
            avrdir    - the average direction computed from neighbors
            dir_strength - the strength of the average direction
            nvalid    - the number of valid directions used to compute the
                        average
      **************************************************************************/
      void Average8NeighDir(int *avrdir, double *dir_strength, int *nvalid,
                           int *imap, const int mx, const int my,
                           const int mw, const int mh,
                           const Dir2Rad& dir2rad);

      /*************************************************************************
      **************************************************************************
      TestTopEdge - Walks the top edge of a concentric square in the IMAP,
                     testing directions along the way to see if they should
                     be removed due to being too weak or inconsistent with
                     respect to their adjacent neighbors.

         Input:
            lbox      - left edge of current concentric square
            tbox      - top edge of current concentric square
            rbox      - right edge of current concentric square
            imap      - vector of IMAP integer directions
            mw        - width (in blocks) of the IMAP
            mh        - height (in blocks) of the IMAP
            dir2rad   - lookup table for converting integer directions
            lqmParams - parameters and thresholds for controlling LQM
         Return Code:
            Positive - direction should be removed from IMAP
            Zero     - direction should NOT be remove from IMAP
      **************************************************************************/
      int TestTopEdge(const int lbox, const int tbox, const int rbox,
                        int *imap, const int mw, const int mh,
                        const Dir2Rad& dir2rad, const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      TestRightEdge - Walks the right edge of a concentric square in the
                     IMAP, testing directions along the way to see if they
                     should be removed due to being too weak or inconsistent
                     with respect to their adjacent neighbors.

         Input:
            tbox      - top edge of current concentric square
            rbox      - right edge of current concentric square
            bbox      - bottom edge of current concentric square
            imap      - vector of IMAP integer directions
            mw        - width (in blocks) of the IMAP
            mh        - height (in blocks) of the IMAP
            dir2rad   - lookup table for converting integer directions
            lqmParams - parameters and thresholds for controlling LQM
         Return Code:
            Positive - direction should be removed from IMAP
            Zero     - direction should NOT be remove from IMAP
      **************************************************************************/
      int TestRightEdge(const int tbox, const int rbox,
                        const int bbox, int *imap, const int mw, const int mh,
                        const Dir2Rad& dir2rad, const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      TestBottomEdge - Walks the bottom edge of a concentric square in the
                  IMAP, testing directions along the way to see if they
                  should be removed due to being too weak or inconsistent
                  with respect to their adjacent neighbors.
         Input:
            lbox      - left edge of current concentric square
            rbox      - right edge of current concentric square
            bbox      - bottom edge of current concentric square
            imap      - vector of IMAP integer directions
            mw        - width (in blocks) of the IMAP
            mh        - height (in blocks) of the IMAP
            dir2rad   - lookup table for converting integer directions
            lqmParams - parameters and thresholds for controlling LQM
         Return Code:
            Positive - direction should be removed from IMAP
            Zero     - direction should NOT be remove from IMAP
      **************************************************************************/
      int TestBottomEdge(const int lbox, const int rbox,
                           const int bbox, int *imap, const int mw, const int mh,
                           const Dir2Rad& dir2rad, const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      TestLeftEdge - Walks the left edge of a concentric square in the IMAP,
                     testing directions along the way to see if they should
                     be removed due to being too weak or inconsistent with
                     respect to their adjacent neighbors.

         Input:
            lbox      - left edge of current concentric square
            tbox      - top edge of current concentric square
            bbox      - bottom edge of current concentric square
            imap      - vector of IMAP integer directions
            mw        - width (in blocks) of the IMAP
            mh        - height (in blocks) of the IMAP
            dir2rad   - lookup table for converting integer directions
            lqmParams - parameters and thresholds for controlling LQM
         Return Code:
            Positive - direction should be removed from IMAP
            Zero     - direction should NOT be remove from IMAP
      **************************************************************************/
      int TestLeftEdge(const int lbox, const int tbox,
                        const int bbox, int *imap, const int mw, const int mh,
                        const Dir2Rad& dir2rad, const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      SmoothDirectionMap - Takes a vector of integer directions and smooths
                  them by analyzing the direction of adjacent neighbors.

         Input:
            direction_map - (Input/Output) vector of integer block values, smoothed in-place
            mw        - width (in blocks) of the map
            mh        - height (in blocks) of the map
            dir2rad   - lookup table for converting integer directions
            lqmParams - parameters and thresholds for controlling LQM
      **************************************************************************/
      void SmoothDirectionMap(int *direction_map, int *low_contrast_map,
                     const int mw, const int mh,
                     const Dir2Rad& dir2rad, const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      InterpolateDirectionMap - Take a Direction Map and Low Contrast
                  Map and attempts to fill in INVALID directions in the
                  Direction Map based on a blocks valid neighbors.  The
                  valid neighboring directions are combined in a weighted
                  average inversely proportional to their distance from
                  the block being interpolated.  Low Contrast blocks are
                  used to prempt the search for a valid neighbor in a
                  specific direction, which keeps the process from
                  interpolating directions for blocks in the background and
                  and perimeter of the fingerprint in the image.

         Input:
            direction_map    - (Input/Output) map of blocks containing directional ridge flow, overwritten by interpolated map
            low_contrast_map - map of blocks flagged as LOW CONTRAST
            mw        - number of blocks horizontally in the maps
            mh        - number of blocks vertically in the maps
            lqmParams - parameters and thresholds for controlling LQM
      **************************************************************************/
      void InterpolateDirectionMap(int *direction_map, int *low_contrast_map,
                           const int mw, const int mh, const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      FindValidBlock - Take a Direction Map, Low Contrast Map,
                  Starting block address, a direction and searches the
                  maps in the specified direction until either a block valid
                  direction is encountered or a block flagged as LOW CONTRAST
                  is encountered.  If a valid direction is located, it and the
                  address of the corresponding block are returned with a
                  code of FOUND.  Otherwise, a code of NOT_FOUND is returned.

         Input:
            direction_map    - map of blocks containing directional ridge flows
            low_contrast_map - map of blocks flagged as LOW CONTRAST
            sx        - X-block coord where search starts in maps
            sy        - Y-block coord where search starts in maps
            mw        - number of blocks horizontally in the maps
            mh        - number of blocks vertically in the maps
            x_incr    - X-block increment to direct search
            y_incr    - Y-block increment to direct search
         Output:
            nbr_dir   - valid direction found
            nbr_x     - X-block coord where valid direction found
            nbr_y     - Y-block coord where valid direction found
         Return Code:
            LQM_FOUND     - neighboring block with valid direction found
            LQM_NOT_FOUND - neighboring block with valid direction NOT found
      **************************************************************************/
      int FindValidBlock(int *nbr_dir, int *nbr_x, int *nbr_y,
                           int *direction_map, int *low_contrast_map,
                           const int sx, const int sy,
                           const int mw, const int mh,
                           const int x_incr, const int y_incr);

      /*************************************************************************
      **************************************************************************
      SetMarginBlocks - Take an image map and sets its perimeter values to
                  the specified value.

         Input:
            map       - map of blocks to be modified
            mw        - number of blocks horizontally in the map
            mh        - number of blocks vertically in the map
            margin_value - value to be assigned to the perimeter blocks
         Output:
            map       - resulting map
      **************************************************************************/
      void SetMarginBlocks(int *map, const int mw, const int mh,
                        const int margin_value);

      /*************************************************************************
      **************************************************************************
      GenHighCurveMap - Takes a Direction Map and generates a new map
               that flags blocks with HIGH CURVATURE.

         Input:
            direction_map - map of blocks containing directional ridge flow
            mw        - the width (in blocks) of the map
            mh        - the height (in blocks) of the map
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            out_hcmap    - points to the created High Curvature Map
      **************************************************************************/
      void GenHighCurveMap(std::vector<int>& out_hcmap, int *direction_map,
                        const int mw, const int mh, const LQMParams& lqmParams,
                        unsigned char *validNeighbors, unsigned char *curvatureMap);

      /*************************************************************************
      **************************************************************************
      NumValid8Neigh - Given a block in an IMAP, counts the number of
                        immediate neighbors that have a valid IMAP direction.

         Input:
            imap - 2-D vector of directional ridge flows
            mx   - horizontal coord of current block in IMAP
            my   - vertical coord of current block in IMAP
            mw   - width (in blocks) of the IMAP
            mh   - height (in blocks) of the IMAP
         Return Code:
            Non-negative - the number of valid IMAP neighbors
      **************************************************************************/
      int NumValid8Neigh(int *imap, const int mx, const int my,
                        const int mw, const int mh);

      /*************************************************************************
      **************************************************************************
      Vorticity - Measures the amount of cumulative curvature incurred
                  among the IMAP neighbors of the given block.

         Input:
            imap  - 2D vector of ridge flow directions
            mx    - horizontal coord of current IMAP block
            my    - vertical coord of current IMAP block
            mw    - width (in blocks) of the IMAP
            mh    - height (in blocks) of the IMAP
            ndirs - number of possible directions in the IMAP
         Return Code:
            Non-negative - the measured vorticity among the neighbors
      **************************************************************************/
      int vorticity(int *imap, const int mx, const int my,
                  const int mw, const int mh, const int ndirs);

      /*************************************************************************
      **************************************************************************
      AccumNeighVorticity - Accumlates the amount of curvature measures
                           between neighboring IMAP blocks.

         Input:
            dir1  - first neighbor's integer IMAP direction
            dir2  - second neighbor's integer IMAP direction
            ndirs - number of possible IMAP directions
         Output:
            vmeasure - accumulated vorticity among neighbors measured so far
      **************************************************************************/
      void AccumNeighVorticity(int *vmeasure, const int dir1, const int dir2,
                              const int ndirs);

      /*************************************************************************
      **************************************************************************
      curvature - Measures the largest change in direction between the
                  current IMAP direction and its immediate neighbors.

         Input:
            imap  - 2D vector of ridge flow directions
            mx    - horizontal coord of current IMAP block
            my    - vertical coord of current IMAP block
            mw    - width (in blocks) of the IMAP
            mh    - height (in blocks) of the IMAP
            ndirs - number of possible directions in the IMAP
         Return Code:
            Non-negative - maximum change in direction found (curvature)
            Negative     - No valid neighbor found to measure change in direction
      **************************************************************************/
      int curvature(int *imap, const int mx, const int my,
                  const int mw, const int mh, const int ndirs);

      /*************************************************************************
      **************************************************************************
      ClosestDirDist - Takes to integer IMAP directions and determines the
                        closest distance between them accounting for
                        wrap-around either at the beginning or ending of
                        the range of directions.

         Input:
            dir1  - integer value of the first direction
            dir2  - integer value of the second direction
            ndirs - the number of possible directions
         Return Code:
            Non-negative - distance between the 2 directions
      **************************************************************************/
      int ClosestDirDist(const int dir1, const int dir2, const int ndirs);

      void GenDirectionChangeMap(unsigned char *directionChangeMap, int *initialDirectionMap, int *direction_map, int mw, int mh);

      /*************************************************************************
      **************************************************************************
      Binarize - Takes a padded grayscale input image and its associated
                  Direction Map and produces a binarized version of the
                  image.  It then fills horizontal and vertical "holes" in
                  the binary image results.  Note that the input image must
                  be padded sufficiently to contain in memory rotated
                  directional binarization grids applied to pixels along the
                  perimeter of the input image.

         Input:
            pdata       - padded input grayscale image
            pw          - padded width (in pixels) of input image
            ph          - padded height (in pixels) of input image
            direction_map - 2-D vector of discrete ridge flow directions
            mw          - width (in blocks) of the map
            dirbingrids - set of rotated grid offsets used for directional
                        binarization
            lqmParams   - parameters and thresholds for controlling LQM
         Output:
            odata - reference to created (unpadded) binary image
            ow    - width of binary image
            oh    - height of binary image
      **************************************************************************/
      void Binarize(std::vector<unsigned char>& odata, int& ow, int& oh,
               unsigned char *pdata, const int pw, const int ph,
               int *direction_map, const int mw,
               const RotGrids& dirbingrids, const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      BinarizeImage - Takes a grayscale input image and its associated
                  Direction Map and generates a binarized version of the
                  image.  Note that there is no "Isotropic" binarization
                  used in this version.

         Input:
            pdata       - padded input grayscale image
            pw          - padded width (in pixels) of input image
            ph          - padded height (in pixels) of input image
            direction_map - 2-D vector of discrete ridge flow directions
            mw          - width (in blocks) of the map
            blocksize   - dimension (in pixels) of each NMAP block
            dirbingrids - set of rotated grid offsets used for directional
                        binarization
         Output:
            odata  - points to binary image results
            ow     - points to binary image width
            oh     - points to binary image height
      **************************************************************************/
      void BinarizeImage(std::vector<unsigned char>& odata, int& ow, int& oh,
                        unsigned char *pdata, const int pw, const int ph,
                        const int *direction_map, const int mw,
                        const int blocksize, const RotGrids& dirbingrids);

      /*************************************************************************
      **************************************************************************
      BinarizeFromDirection - Determines the binary value of a grayscale pixel based
                  on a VALID IMAP ridge flow direction.

         CAUTION: The image to which the input pixel points must be appropriately
                  padded to account for the radius of the rotated grid.  Otherwise,
                  this routine may access "unkown" memory.

         Input:
            pptr        - pointer to current grayscale pixel
            idir        - IMAP integer direction associated with the block the
                        current is in
            dirbingrids - set of precomputed rotated grid offsets
         Return Code:
            LQM_BLACK_PIXEL - pixel intensity for BLACK
            LQM_WHITE_PIXEL - pixel intensity of WHITE
      **************************************************************************/
      int BinarizeFromDirection(const unsigned char *pptr, const int idir,
                     const RotGrids& dirbingrids);

      /*************************************************************************
      **************************************************************************
      FillHoles - Takes an input image and analyzes triplets of horizontal
                  pixels first and then triplets of vertical pixels, filling
                  in holes of width 1.  A hole is defined as the case where
                  the neighboring 2 pixels are equal, AND the center pixel
                  is different.  Each hole is filled with the value of its
                  immediate neighbors. This routine modifies the input image.

         Input:
            bdata - binary image data to be processed
            iw    - width (in pixels) of the binary input image
            ih    - height (in pixels) of the binary input image
         Output:
            bdata - points to the results
      **************************************************************************/
      void FillHoles(unsigned char *bdata, const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      GrayToBin - Takes an 8-bit threshold value and two 8-bit pixel values.
               Those pixels in the image less than the threhsold are set
               to the first specified pixel value, whereas those pixels
               greater than or equal to the threshold are set to the second
               specified pixel value.  On application for this routine is
               to convert binary images from 8-bit pixels valued {0,255} to
               {1,0} and vice versa.

         Input:
            thresh      - 8-bit pixel threshold
            less_pix    - pixel value used when image pixel is < threshold
            greater_pix - pixel value used when image pixel is >= threshold
            bdata       - 8-bit image data
            iw          - width (in pixels) of the image
            ih          - height (in pixels) of the image
         Output:
            bdata       - altered 8-bit image data
      **************************************************************************/
      void GrayToBin(const int thresh, const int less_pix, const int greater_pix,
                  unsigned char *bdata, const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      PixelizeMap - Takes a block image map and assigns each pixel in the
               image its corresponding block value.  This allows block
               values in maps to be directly accessed via pixel addresses.

         Input:
            iw        - the width (in pixels) of the corresponding image
            ih        - the height (in pixels) of the corresponding image
            imap      - input block image map
            mw        - the width (in blocks) of the map
            mh        - the height (in blocks) of the map
            blocksize - the dimension (in pixels) of each block
         Output:
            omap      - points to the resulting pixelized map
      **************************************************************************/
      void PixelizeMap(std::vector<int>& omap, const int iw, const int ih, int *imap, const int mw, const int mh, const int blocksize);

      /*************************************************************************
      **************************************************************************
      FreePath - Traverses a straight line between 2 pixel points in an
                  image and determines if a "free path" exists between the
                  2 points by counting the number of pixel value transitions
                  between adjacent pixels along the trajectory.

         Input:
            x1       - x-pixel coord of first point
            y1       - y-pixel coord of first point
            x2       - x-pixel coord of second point
            y2       - y-pixel coord of second point
            bdata    - binary image data (0==while & 1==black)
            iw       - width (in pixels) of image
            lqmParams- parameters and thresholds specified by LQM
         Return Code:
            true      - free path determined to exist
            false     - free path determined not to exist
      **************************************************************************/
      bool FreePath(const int x1, const int y1, const int x2, const int y2,
                  unsigned char *bdata, const int iw,
                  const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      LinePoints - Returns the contiguous coordinates of a line connecting
                  2 specified points.

         Input:
            x1      - x-coord of first point
            y1      - y-coord of first point
            x2      - x-coord of second point
            y2      - y-coord of second point
         Output:
            points - list of points along line trajectory
      **************************************************************************/
      void LinePoints(std::vector<std::pair<int, int>>& points, const int x1, const int y1, const int x2, const int y2);

      OpenLQM::Fingerprint CreateNormFingerprint(const OpenLQM::Fingerprint& img);

      template <typename T>
      T IntDiv(T num, T den, bool roundUp) {
         if (roundUp) {
            return (num + den - 1) / den;
         } else {
            return num / den;
         }
      }

      // Round up to the next wordSize
      template <typename T>
      T WordAlign(T val, bool roundUp, T wordSize = 4) {
         return IntDiv(val, wordSize, roundUp) * wordSize;
      }

      float ClampResolution(float ppi);
   }
}
