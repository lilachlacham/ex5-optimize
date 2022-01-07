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
    //int num;
} pixel_sum;


/* Compute min and max of two integers, respectively */
//change min, max functions to MACRO
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))
#define calcIndex(i, j, n) ((i)*(n)+(j))


/*
 * initialize_pixel_sum - Initializes all fields of sum to 0
 */
#define initialize_pixel_sum(sum) (sum)->red = (sum)->green = (sum)->blue = 0

/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 */

#define assign_sum_to_pixel(current_pixel,sum, kernelScale) \
	sum.red = sum.red / kernelScale; \
	sum.green = sum.green / kernelScale; \
	sum.blue = sum.blue / kernelScale; \
	current_pixel.red = (unsigned char) (min(max(sum.red, 0), 255));\
	current_pixel.green = (unsigned char) (min(max(sum.green, 0), 255));\
	current_pixel.blue = (unsigned char) (min(max(sum.blue, 0), 255));


#define filter(ii, jj, loop_pixel) \
 intensity = (int) loop_pixel.red + ((int) loop_pixel.green) + ((int) loop_pixel.blue);\
if (intensity <= min_intensity) {   \
    min_intensity = intensity;  \
    min_row = ii;   \
    min_col = jj;   \
}   \
if (intensity > max_intensity) {    \
    max_intensity = intensity;  \
    max_row = ii;   \
    max_col = jj;   \
}

/*
* sum_pixels_by_weight - Sums pixel values, scaled by given weight
*/
#define sum_pixels_by_weight(sum, p, weight) \
(sum).red += ((int) (p).red) * (weight); \
(sum).green += ((int) (p).green) * (weight);\
(sum).blue += ((int) (p).blue) * (weight);

#define sum_pixels(sum, p) \
(sum).red += ((int) (p).red); \
(sum).green += ((int) (p).green);\
(sum).blue += ((int) (p).blue);

#define sum_pixels_by_color_with_kernel(color, index) \
sum.color = (src[index].color * (-1)) + (src[index+1].color * (-1)) + \
(src[index+2].color * (-1)) + (src[index+dim].color * (-1)) + \
(src[index+dim+1].color * (9)) + (src[index+dim+2].color * (-1)) + \
(src[index+dim+dim].color * (-1)) + (src[index+dim+dim+1].color * (-1)) + \
(src[index+dim+dim+2].color * (-1));

#define sum_pixels_by_color_without_kernel(color, index) \
sum.color = (src[index].color) + (src[index+1].color) + \
(src[index+2].color) + (src[index+dim].color) + \
(src[index+dim+1].color) + (src[index+dim+2].color) + \
(src[index+dim+dim].color) + (src[index+dim+dim+1].color) + \
(src[index+dim+dim+2].color);


/*
* [1, 1, 1]
* [1, 1, 1]
* [1, 1, 1]
*/
//int blurKernel[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

/*
* [-1, -1, -1]
* [-1, 9, -1]
* [-1, -1, -1]
*/
//int sharpKernel[3][3] = {{-1,-1,-1},{-1,9,-1},{-1,-1,-1}};
/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
void smoothSharp(int dim, pixel *src, pixel *dst, int kernelScale, bool filter) {



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
            register int firstSrcIndex = calcIndex(ii, jj, dim);

            sum_pixels_by_color_with_kernel(red, firstSrcIndex);
            sum_pixels_by_color_with_kernel(green, firstSrcIndex);
            sum_pixels_by_color_with_kernel(blue, firstSrcIndex);

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
            assign_sum_to_pixel(current_pixel, sum, kernelScale);
            dst[calcIndex(i, j, dim)] = current_pixel;
        }
    }
}
void smoothBlur(int dim, pixel *src, pixel *dst, int kernelScale, bool filter) {

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
            register int firstSrcIndex = calcIndex(ii, jj, dim);
            sum_pixels_by_color_without_kernel(red, firstSrcIndex);
            sum_pixels_by_color_without_kernel(green, firstSrcIndex);
            sum_pixels_by_color_without_kernel(blue, firstSrcIndex);


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
            assign_sum_to_pixel(current_pixel, sum, kernelScale);
            dst[calcIndex(i, j, dim)] = current_pixel;
        }
    }
}


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

void pixelsToChars(pixel* pixels, pixel* pixel1) {
	register int row, col;
	for (row = 0 ; row < m ; ++row) {
		for (col = 0 ; col < n ; ++col) {
            register unsigned int row_n_col = row*n + col;
            register unsigned int row_n_col_3 = 3*row_n_col;
            *(pixel *)&(image->data[row_n_col_3]) = pixels[row_n_col];
            pixel1[row_n_col]= pixels[row_n_col];
		}
	}
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

    pixel* pixelsImg = malloc(m*n*sizeof(pixel));
    pixel* backupOrg = malloc(m*n*sizeof(pixel));
    charsToPixels2(pixelsImg, backupOrg);

	if (flag == '1') {	
		// blur image
        smoothBlur(m, backupOrg, pixelsImg, 9, false);
        pixelsToChars(pixelsImg, backupOrg);
		// write result image to file
		writeBMP(image, srcImgpName, blurRsltImgName);	

        // sharpen the resulting image
        smoothSharp(m, backupOrg, pixelsImg, 1, false);
        pixelsToChars(pixelsImg, backupOrg);
		// write result image to file
		writeBMP(image, srcImgpName, sharpRsltImgName);	
	} else {
		// apply extermum filtered kernel to blur image
        smoothBlur(m, backupOrg, pixelsImg,7,true);
        pixelsToChars(pixelsImg, backupOrg);

		// write result image to file
		writeBMP(image, srcImgpName, filteredBlurRsltImgName);

		// sharpen the resulting image
        smoothSharp(m,backupOrg,pixelsImg,1,false);
        pixelsToChars(pixelsImg, backupOrg);
		// write result image to file
		writeBMP(image, srcImgpName, filteredSharpRsltImgName);

	}
    free(pixelsImg);
    free(backupOrg);
}

