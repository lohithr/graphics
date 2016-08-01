#include <algorithm>
#include <cassert>
#include <cctype>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <string>
#include <math.h>
#include "stb_image_write.hpp"

static int const WIDTH         =  800;  ///< Width of output image in pixels.
static int const HEIGHT        =  800;  ///< Height of output image in pixels.
static int const NUM_CHANNELS  =    3;  ///< Number of channels (color components) per output pixel.

typedef std::complex<double>  Complex;  ///< A double-precision complex number.
typedef unsigned char         uint8;    ///< A single byte.

/** A representation of a color in the RGB model, with one byte per channel. */
struct RGB8
{
  uint8 r;  ///< Red component, range 0 to 255.
  uint8 g;  ///< Green component, range 0 to 255.
  uint8 b;  ///< Blue component, range 0 to 255.
};

/** Check if the \a pattern string is a suffix of the \a test string. */
bool endsWith(std::string const & test, std::string const & pattern);

/** Return the lowercase version of a string. */
std::string toLower(std::string const & s);

/**
 * Convert pixel coordinates in an image to the corresponding point location in the complex plane.
 *
 * @param row Row containing pixel (in range 0 to HEIGHT - 1).
 * @param col Column containing pixel (in range 0 to WIDTH - 1).
 * @param lo Lower left corner of the box.
 * @param real_range Width (along real axis) of corresponding region in the complex plane.
 * @param imag_range Height (along imaginary axis) of corresponding region in the complex plane.
 *
 * @return The complex number corresponding to the pixel location (row, col).
 */
Complex
pixelToPoint(int row, int col, Complex const & lo, double real_range, double imag_range)
{
  assert(real_range > 0 && imag_range > 0);

  double imag_step = imag_range / HEIGHT, real_step = real_range / WIDTH;
  double _real , _imag;
  _imag = lo.imag() + (imag_step * (HEIGHT - row)) - (0.5 * imag_step); //mapping imaginary part to row in pixel
  _real = lo.real() + (col * real_step) + (0.5 * real_step);  //mapping real part to column in pixel
  return Complex( _real, _imag);
}

/**
 * Check for convergence of the iterated Mandelbrot function and return an appropriate color for the point at complex plane
 * location c = x + iy.
 */
RGB8
getMandelbrotPointColor(Complex const & c)
{
  RGB8 color;
  Complex com=c;
  int _count;
  for ( _count = 0; !(com.imag() > 2.0 || com.imag() < -2.0 || com.real() > 2.0 || com.real() < -2.0); _count++) // exits if point crosses the -2,-2 to 2,2 square boundary
  {
   if(_count==400)
   { 
     //colouring scheme for points inside the mandelbrot set
     color.r = 0;
     color.g = 0;
     color.b = 0;
     return color;
   }
    com=pow(com,2)+c; //applying the funtion f(com) = com^2 + c 
  } 
  //colouring scheme for points outside the mandelbrot set
  color.r =  ( 102 * _count ) / 40 ;  
  color.g = ( 255 * _count ) / 40 ;
  color.b = (153 * _count) / 40 ; 

  return color;
}

/**
 * Generates an image of the Mandelbrot set in the complex plane, within a square box.
 *
 * @param lo Lower left corner of the box.
 * @param s Length of a side of the square box.
 * @param out_path The path to the file in which to write the output (should end with .png).
 *
 * @return True on success, false on error.
 */
bool
generateMandelbrotImage(Complex const & lo, double s, std::string const & out_path)
{
  // Sanity checks
  if (s <= 0)
  {
    std::cerr << "Side length must be positive" << std::endl;
    return false;
  }

  if (out_path.empty() || !endsWith(toLower(out_path), ".png"))
  {
    std::cerr << "Output format must be PNG" << std::endl;
    return false;
  }

  // Allocate a buffer for the image data. Pixels are stored row-by-row: left to right, top to bottom. For example, an image
  // looking like this:
  //
  // o o o o o
  // o * o o o
  // o * o o o
  // o * * * o
  // o o o o o
  //
  // Would be stored as:
  //
  // o o o o o o * o o o o * o o o o * * * o o o o o o
  //
  // Each pixel has three byte-sized components: red, green and blue, in that order.
  //
  long num_bytes = WIDTH * HEIGHT * NUM_CHANNELS;
  uint8 * data = (uint8 *)std::malloc((std::size_t)num_bytes);  // assumes stbi uses malloc
  if (!data)
  {
    std::cerr << "Could not allocate image buffer" << std::endl;
    return false;
  }

  // Loop over the pixels, computing the color of each one based on whether it is in the Mandelbrot set or not.
  for (int row = 0; row < HEIGHT; row++)
  {
    for (int col = 0; col < WIDTH; col++)
    {
      Complex c = pixelToPoint(row, col, lo, s, s);
      RGB8 color = getMandelbrotPointColor(c);

      long index = (row * WIDTH + col) * NUM_CHANNELS;
      data[index    ] = color.r;
      data[index + 1] = color.g;
      data[index + 2] = color.b;
    }
  }

  // Write image buffer to the output file
  if (!stbi_write_png(out_path.c_str(), WIDTH, HEIGHT, NUM_CHANNELS, data, WIDTH * NUM_CHANNELS))
  {
    std::cerr << "Could not write output image" << std::endl;
    return false;
  }

  return true;
}

int
main(int argc, char * argv[])
{
  if (argc != 5)
  {
    std::cerr << "Usage: " << argv[0] << " <lo.real>  <lo.imag>  <s>  <out_path>" << std::endl;
    return -1;
  }

  double lo_real        =  std::atof(argv[1]);
  double lo_imag        =  std::atof(argv[2]);
  double s              =  std::atof(argv[3]);
  std::string out_path  =  argv[4];

  if (!generateMandelbrotImage(Complex(lo_real, lo_imag), s, out_path))
    return -1;

  return 0;
}

//============================================================================================================================
// Implementations of utility functions
//============================================================================================================================

bool
endsWith(std::string const & test, std::string const & pattern)
{
  if (test.size() >= pattern.size())
  {
    long te = (long)test.size() - 1;
    long pe = (long)pattern.size() - 1;

    for (long i = (long)pattern.size() - 1; i >= 0; --i)
    {
      if (pattern[(size_t)(pe - i)] != test[(size_t)(te - i)])
        return false;
    }

    return true;
  }
  else
    return false;
}

std::string
toLower(std::string const & s)
{
  std::string result = s;
  std::transform(result.begin(), result.end(), result.begin(), ::tolower);
  return result;
}
