//207375700 Racheli Lilach Lacham
#include <stdbool.h>

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
//change min, max , calcIndex functions to MACRO
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))
#define calcIndex(i, j, n) ((i)*(n)+(j))


/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 * change the function to macro
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
(sum).red += ((int) (p)->red) * (weight); \
(sum).green += ((int) (p)->green) * (weight);\
(sum).blue += ((int) (p)->blue) * (weight);


#define sum_pixels_by_color_with_kernel(color, index) \
sum.color = (src[index].color * (-1)) + (src[index+1].color * (-1)) + \
(src[index+2].color * (-1)) + (src[index+dim].color * (-1)) + \
(src[index+dim+1].color * (9)) + (src[index+dim+2].color * (-1)) + \
(src[index+dim+dim].color * (-1)) + (src[index+dim+dim+1].color * (-1)) + \
(src[index+dim+dim+2].color * (-1));

#define sum_pixels_by_color_without_kernel(sum, color, index) \
sum.color += (src[index].color) + (src[(index)+dim].color) + \
(src[(index)+dim+dim].color);

#define sumPixel(ptrRow, weight) \
sum.red += ptrRow->red * weight; \
sum.green += ptrRow->green * weight; \
sum.blue += ptrRow->blue * weight; \
if (filter) { \
    intensity = (int) ptrRow->red + ((int) ptrRow->green) + ((int) ptrRow->blue); \
    if (intensity <= min_intensity) { \
    min_intensity = intensity; \
    minPtr = ptrRow; \
} \
    if (intensity > max_intensity) { \
    max_intensity = intensity; \
    maxPtr = ptrRow; \
} \
} \
++ptrRow;

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
/*
 * create one function for sharp and another for blur, so no need to get as parameter the kernel.
 */
void smoothSharp(int dim, pixel *src, pixel *dst, int kernelScale, bool filter) {
    register int i;
    register int end = dim - 1;
    register int intensity;
    pixel* minPtr;
    pixel* maxPtr;
    pixel* ptrRow1 = src;
    pixel* ptrRow2 = src + dim;
    pixel* ptrRow3 = src + dim + dim;
    dst += dim+1;

    for (i = 1 ; i < end; ++i) {
        register int j = 1;
        pixel* ptrEndRow1 = ptrRow1 + dim;
        while(ptrRow1 < ptrEndRow1) {
            pixel_sum sum = {0,0,0};
            pixel current_pixel;
            register int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
            register int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0

            sumPixel(ptrRow1, -1);
            sumPixel(ptrRow1, -1);
            sumPixel(ptrRow1, -1);
            ptrRow1 -=2;
            sumPixel(ptrRow2, -1);
            sumPixel(ptrRow2, 9);
            sumPixel(ptrRow2, -1);
            ptrRow2 -=2;
            sumPixel(ptrRow3, -1);
            sumPixel(ptrRow3, -1);
            sumPixel(ptrRow3, -1);
            ptrRow3 -=2;

            if (filter) {
                // find min and max coordinates
                // filter out min and max
                sum_pixels_by_weight(sum,minPtr, -1);
                sum_pixels_by_weight(sum, maxPtr, -1);
            }
            assign_sum_to_pixel(current_pixel, sum, kernelScale);
            *dst = current_pixel;
            ++dst;
        }
    }
}

void smoothBlur(int dim, pixel *src, pixel *dst, int kernelScale, bool filter) {
    register int i, j;
    register int end = dim - 1;
    register int firstSrcIndex;
    for (i = 1 ; i < end; ++i) {
        pixel_sum sum1 = {0,0,0};
        pixel_sum sum2 = {0,0,0};
        int ii = max(i-1, 0);
        firstSrcIndex = calcIndex(ii, 0, dim);
        //calculate for every new mask the two left Columns.
        sum_pixels_by_color_without_kernel(sum1, red, firstSrcIndex);
        sum_pixels_by_color_without_kernel(sum2, red, firstSrcIndex+1);
        sum_pixels_by_color_without_kernel(sum1, green, firstSrcIndex);
        sum_pixels_by_color_without_kernel(sum2, green, firstSrcIndex+1);
        sum_pixels_by_color_without_kernel(sum1, blue, firstSrcIndex);
        sum_pixels_by_color_without_kernel(sum2, blue, firstSrcIndex+1);
        for (j = 1 ; j < end ; ++j, ++firstSrcIndex) {
            pixel_sum sum = {0,0,0};
            pixel current_pixel;
            register int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
            register int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
            register int min_row, min_col, max_row, max_col;
            pixel loop_pixel;

            register int jj = max(j-1, 0);

            pixel_sum tempSum = {0,0,0};
            //calculate the right Column and saved it for the next loop.
            sum_pixels_by_color_without_kernel(sum, red, firstSrcIndex+2);
            sum_pixels_by_color_without_kernel(sum, green, firstSrcIndex+2);
            sum_pixels_by_color_without_kernel(sum, blue, firstSrcIndex+2);
            tempSum = sum;
            sum.red += sum1.red + sum2.red;
            sum.green += sum1.green + sum2.green;
            sum.blue += sum1.blue + sum2.blue;
            sum1 = sum2;
            sum2 = tempSum;

            if (filter) {
                register int maxII = min(i+1, dim-1);
                register int maxJJ = min(j+1, dim-1);
                for(ii = max(i-1, 0); ii <= maxII; ++ii) {
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
                sum_pixels_by_weight(sum, &(src[calcIndex(min_row, min_col, dim)]), -1);
                sum_pixels_by_weight(sum, &(src[calcIndex(max_row, max_col, dim)]), -1);
            }

            // assign kernel's result to pixel at [i,j]
            assign_sum_to_pixel(current_pixel, sum, kernelScale);
            dst[calcIndex(i, j, dim)] = current_pixel;
        }
    }
}


#define charsToPixels2(pixels, pixels2) \
    /*do{                                    \
    register int row, col;  */            \
    memcpy(pixels,image->data,m*n*sizeof(pixel));                                   \
    memcpy(pixels2,pixels,m*n*sizeof(pixel));                                   \
    /*for (row = m ; row-- ;) { \
        for (col = n ; col--;) { \
            register unsigned int row_n_col = row*n + col; \
            register unsigned int row_n_col_3 = (row_n_col<<1) + row_n_col; \
            pixels[row_n_col] = *(pixel*)&(image->data[row_n_col_3]); \
            pixels2[row_n_col] = *(pixel*)&(image->data[row_n_col_3]); \
        } \
    }  \
} while(0)*/

#define pixelsToChars(pixels, pixel1) \
    memcpy(image->data,pixels,m*n*sizeof(pixel));                                   \
    memcpy(pixel1,pixels,m*n*sizeof(pixel));                                  \
    /*do{                   \
	register int row, col; \
	for (row = m ; row-- ;) { \
		for (col = n ; col-- ;) { \
            register unsigned int row_n_col = row*n + col; \
            register unsigned int row_n_col_3 = (row_n_col<<1) + row_n_col; \
            *(pixel *)&(image->data[row_n_col_3]) = pixels[row_n_col]; \
            pixel1[row_n_col]= pixels[row_n_col]; \
		} \
	} \
} while(0)*/

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

