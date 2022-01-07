#include <stdbool.h> 
#define KERNEL_SIZE 3


typedef struct {
   unsigned char red;
   unsigned char green;
   unsigned char blue;
} pixel;

typedef struct {
    int red;
    int green;
    int blue;
    // int num;
} pixel_sum;


/* Compute min and max of two integers, respectively */
//change min, max functions to MACRO
#define min(a, b) (a < b ? a : b)
#define max(a, b) (a > b ? a : b)
#define calcIndex(i, j, n) ((i)*(n)+(j))

#define filterMacro(src, ii, jj, dim, min_intensity, min_row, min_col, max_intensity, max_row, max_col) \
loop_pixel = src[calcIndex(ii, jj, dim)];\
intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));\
if (intensity <= min_intensity) {\
    min_intensity = intensity;\
    min_row = ii;\
    min_col = jj;\
}\
if (intensity > max_intensity) {\
    max_intensity = intensity;\
    max_row = ii;\
    max_col = jj;\
}

/*
 * initialize_pixel_sum - Initializes all fields of sum to 0
 */
#define initialize_pixel_sum(sum) (sum)->red = (sum)->green = (sum)->blue = 0

/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) {

	// divide by kernel's weight
	sum.red = sum.red / kernelScale;
	sum.green = sum.green / kernelScale;
	sum.blue = sum.blue / kernelScale;

	// truncate each pixel's color values to match the range [0,255]
	current_pixel->red = (unsigned char) (min(max(sum.red, 0), 255));
	current_pixel->green = (unsigned char) (min(max(sum.green, 0), 255));
	current_pixel->blue = (unsigned char) (min(max(sum.blue, 0), 255));
	return;
}

/*
* sum_pixels_by_weight - Sums pixel values, scaled by given weight
*/
#define sum_pixels_by_weight(sum, p, weight) \
(sum).red += ((int) (p).red) * (weight); \
(sum).green += ((int) (p).green) * (weight);\
(sum).blue += ((int) (p).blue) * (weight);

#define sum_pixels_by_color(color) \
sum.color = (src[calcIndex(ii, jj, dim)].color * kernel1) + (src[calcIndex(ii, jj+1, dim)].color * kernel2) + \
(src[calcIndex(ii, jj+2, dim)].color * kernel3) + (src[calcIndex(ii+1, jj, dim)].color * kernel4) + \
(src[calcIndex(ii+1, jj+1, dim)].color * kernel5) + (src[calcIndex(ii+1, jj+2, dim)].color * kernel6) + \
(src[calcIndex(ii+2, jj, dim)].color * kernel7) + (src[calcIndex(ii+2, jj+1, dim)].color * kernel8) + \
(src[calcIndex(ii+2, jj+2, dim)].color * kernel9);

/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
/*void smooth(int dim, pixel *src, pixel *dst, int kernel[KERNEL_SIZE][KERNEL_SIZE], int kernelScale, bool filter) {

    register int kernel1 = kernel[0][0];
    register int kernel2 = kernel[0][1];
    register int kernel3 = kernel[0][2];
    register int kernel4 = kernel[1][0];
    register int kernel5 = kernel[1][1];
    register int kernel6 = kernel[1][2];
    register int kernel7 = kernel[2][0];
    register int kernel8 = kernel[2][1];
    register int kernel9 = kernel[2][2];

    register int i, j;
    register int startBound = KERNEL_SIZE / 2;
    register int endBoundI = dim - startBound;
    register int carry = endBoundI % 20;
    register int endBoundJ = endBoundI - carry;

	for (i = startBound ; i < endBoundI; ++i) {
		for (j =  startBound ; j < endBoundJ ; j+=20) {
#define smoothLoop(j1) \
            do {          \
            pixel_sum sum;\
            pixel current_pixel;\
            register int min_intensity = 766; *//* arbitrary value that is higher than maximum possible intensity, which is 255*3=765*//* \
            register int max_intensity = -1; *//* arbitrary value that is lower than minimum possible intensity, which is 0*//* \
            register int min_row, min_col, max_row, max_col; \
            pixel loop_pixel; \
            initialize_pixel_sum(&sum);\
            register int ii = max(i-1, 0);\
            register int jj = max(j1-1, 0);\
            sum_pixels_by_weight(&sum, src[calcIndex(ii, jj, dim)], kernel1); \
            sum_pixels_by_weight(&sum, src[calcIndex(ii, jj+1, dim)], kernel2); \
            sum_pixels_by_weight(&sum, src[calcIndex(ii, jj+2, dim)], kernel3); \
            sum_pixels_by_weight(&sum, src[calcIndex(ii+1, jj, dim)], kernel4); \
            sum_pixels_by_weight(&sum, src[calcIndex(ii+1, jj+1, dim)], kernel5); \
            sum_pixels_by_weight(&sum, src[calcIndex(ii+1, jj+2, dim)], kernel6); \
            sum_pixels_by_weight(&sum, src[calcIndex(ii+2, jj, dim)], kernel7); \
            sum_pixels_by_weight(&sum, src[calcIndex(ii+2, jj+1, dim)], kernel8); \
            sum_pixels_by_weight(&sum, src[calcIndex(ii+2, jj+2, dim)], kernel9); \
            if (filter) { \
                *//* find min and max coordinates*//* \
                ii = max(i-1, 0); \
                jj = max(j1-1,0); \
                register int intensity = 0; \
                filterMacro(src, ii, jj,dim,min_intensity, min_row, min_col, max_intensity, max_row, max_col) \
                filterMacro(src, ii, jj+1,dim,min_intensity, min_row, min_col, max_intensity, max_row, max_col) \
                filterMacro(src, ii, jj+2,dim,min_intensity, min_row, min_col, max_intensity, max_row, max_col) \
                filterMacro(src, ii+1, jj,dim,min_intensity, min_row, min_col, max_intensity, max_row, max_col) \
                filterMacro(src, ii+1, jj+1,dim,min_intensity, min_row, min_col, max_intensity, max_row, max_col) \
                filterMacro(src, ii+1,jj+2,dim,min_intensity, min_row, min_col, max_intensity, max_row, max_col) \
                filterMacro(src, ii+2, jj,dim,min_intensity, min_row, min_col, max_intensity, max_row, max_col) \
                filterMacro(src, ii+2, jj+1,dim,min_intensity, min_row, min_col, max_intensity, max_row, max_col) \
                filterMacro(src, ii+2,jj+2,dim,min_intensity, min_row, min_col, max_intensity, max_row, max_col) \
                *//* filter out min and max*//* \
                sum_pixels_by_weight(&sum, src[calcIndex(min_row, min_col, dim)], -1);\
                sum_pixels_by_weight(&sum, src[calcIndex(max_row, max_col, dim)], -1);\
            }\
            *//* assign kernel's result to pixel at [i,j]*//* \
            assign_sum_to_pixel(&current_pixel, sum, kernelScale);\
            dst[calcIndex(i, j1, dim)] = current_pixel;\
        } while(0);

            smoothLoop(j);
            smoothLoop(j+1);
            smoothLoop(j+2);
            smoothLoop(j+3);
            smoothLoop(j+4);
            smoothLoop(j+5);
            smoothLoop(j+6);
            smoothLoop(j+7);
            smoothLoop(j+8);
            smoothLoop(j+9);
            smoothLoop(j+10);
            smoothLoop(j+11);
            smoothLoop(j+12);
            smoothLoop(j+13);
            smoothLoop(j+14);
            smoothLoop(j+15);
            smoothLoop(j+16);
            smoothLoop(j+17);
            smoothLoop(j+18);
            smoothLoop(j+19);
		}
        for (;j < endBoundI;++j) {
            smoothLoop(j);
        }
	}
}*/

/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
void smooth(int dim, pixel *src, pixel *dst, int kernel[KERNEL_SIZE][KERNEL_SIZE], int kernelScale, bool filter) {

    register int kernel1 = kernel[0][0];
    register int kernel2 = kernel[0][1];
    register int kernel3 = kernel[0][2];
    register int kernel4 = kernel[1][0];
    register int kernel5 = kernel[1][1];
    register int kernel6 = kernel[1][2];
    register int kernel7 = kernel[2][0];
    register int kernel8 = kernel[2][1];
    register int kernel9 = kernel[2][2];

    register int i, j;
    register int start = KERNEL_SIZE / 2;
    register int end = dim - start;
    for (i = start ; i < end; i++) {
        for (j =  start ; j < end ; j++) {
            pixel_sum sum;
            pixel current_pixel;
            register int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
            register int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
            register int min_row, min_col, max_row, max_col;
            pixel loop_pixel;

            initialize_pixel_sum(&sum);
            register int ii = max(i-1, 0);
            register int jj = max(j-1, 0);

            sum_pixels_by_color(red);
            sum_pixels_by_color(green);
            sum_pixels_by_color(blue);

            if (filter) {
                register int maxII = min(i+1, dim-1);
                register int maxJJ = min(j+1, dim-1);
                for(; ii <= maxII; ++ii) {
                    for(jj = max(j-1, 0); jj <=maxJJ; ++jj) {
                        // apply kernel on pixel at [ii,jj]
                        loop_pixel = src[calcIndex(ii, jj, dim)];
                        register int intensity = (int) loop_pixel.red + ((int) loop_pixel.green) + ((int) loop_pixel.blue);
                        if (intensity <= min_intensity) {
                            min_intensity = intensity;
                            min_row = ii;
                            min_col = jj;
                        }
                        if (intensity > max_intensity) {
                            max_intensity = intensity;
                            max_row = ii;
                            max_col = jj;
                        }
                    }
                }
                // find min and max coordinates
                // filter out min and max
                sum_pixels_by_weight(sum, src[calcIndex(min_row, min_col, dim)], -1);
                sum_pixels_by_weight(sum, src[calcIndex(max_row, max_col, dim)], -1);
            }

            // assign kernel's result to pixel at [i,j]
            assign_sum_to_pixel(&current_pixel, sum, kernelScale);
            dst[calcIndex(i, j, dim)] = current_pixel;
        }
    }
}

/*void charsToPixels(pixel* pixels) {

	register int row, col;
	for (row = 0 ; row < m ; ++row) {
		for (col = 0 ; col < n ; ++col) {
            register unsigned int row_n_col = row*n + col;
            register unsigned int row_n_col_3 = 3*row_n_col;
            pixels[row_n_col] = *(pixel*)&(image->data[row_n_col_3]);
		}
	}
}*/

void charsToPixels2(pixel* pixels, pixel* pixels2) {
    register int row, col;
    for (row = 0 ; row < m ; ++row) {
        for (col = 0 ; col < n ; ++col) {
            register unsigned int row_n_col = row*n + col;
            register unsigned int row_n_col_3 = 3*row_n_col;
            pixels[row_n_col] = *(pixel*)&(image->data[row_n_col_3]);
            pixels2[row_n_col] = *(pixel*)&(image->data[row_n_col_3]);
        }
    }
}

void pixelsToChars(pixel* pixels) {
	register int row, col;
	for (row = 0 ; row < m ; ++row) {
		for (col = 0 ; col < n ; ++col) {
            register unsigned int row_n_col = row*n + col;
            register unsigned int row_n_col_3 = 3*row_n_col;
            *(pixel *)&(image->data[row_n_col_3]) = pixels[row_n_col];
		}
	}
}

/*void copyPixels(pixel* src, pixel* dst) {
	register int row, col;
	for (row = 0 ; row < m ; ++row) {
		for (col = 0 ; col < n ; ++col) {
            register unsigned int row_n_col = row*n + col;
            dst[row_n_col] = src[row_n_col];
		}
	}
}*/

void doConvolution(int kernel[KERNEL_SIZE][KERNEL_SIZE], int kernelScale, bool filter) {

	pixel* pixelsImg = malloc(m*n*sizeof(pixel));
	pixel* backupOrg = malloc(m*n*sizeof(pixel));

	charsToPixels2(pixelsImg, backupOrg);
	//copyPixels(pixelsImg, backupOrg);

	smooth(m, backupOrg, pixelsImg, kernel, kernelScale, filter);

	pixelsToChars(pixelsImg);

	free(pixelsImg);
	free(backupOrg);
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

	/*
	* [1, 1, 1]
	* [1, 1, 1]
	* [1, 1, 1]
	*/
	int blurKernel[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

	/*
	* [-1, -1, -1]
	* [-1, 9, -1]
	* [-1, -1, -1]
	*/
	int sharpKernel[3][3] = {{-1,-1,-1},{-1,9,-1},{-1,-1,-1}};

	if (flag == '1') {	
		// blur image
		doConvolution(blurKernel, 9, false);

		// write result image to file
		writeBMP(image, srcImgpName, blurRsltImgName);	

		// sharpen the resulting image
		doConvolution(sharpKernel, 1, false);
		
		// write result image to file
		writeBMP(image, srcImgpName, sharpRsltImgName);	
	} else {
		// apply extermum filtered kernel to blur image
		doConvolution(blurKernel, 7, true);

		// write result image to file
		writeBMP(image, srcImgpName, filteredBlurRsltImgName);

		// sharpen the resulting image
		doConvolution(sharpKernel, 1, false);

		// write result image to file
		writeBMP(image, srcImgpName, filteredSharpRsltImgName);	
	}
}

