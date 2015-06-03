/* gopt-usage.c version 8.1: tom.viza@gmail.com PUBLIC DOMAIN 2003-8 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "gopt.h"

int main( int argc, const char **argv ){

  FILE *out_stream;
  const char *filename;
  int verbosity;
  int i;  
  void *options= gopt_sort( & argc, argv, gopt_start(
      gopt_option( 'h', 0, gopt_shorts( 'h', '?' ), gopt_longs( "help", "HELP" )),
      gopt_option( 'z', 0, gopt_shorts( 0 ), gopt_longs( "version" )),
      gopt_option( 'v', GOPT_REPEAT, gopt_shorts( 'v' ), gopt_longs( "verbose" )),
      gopt_option( 'o', GOPT_ARG, gopt_shorts( 'o' ), gopt_longs( "output" ))));
  /*
   * there are four possible options to this program, some of which have multiple names:
   *
   * -h -? --help --HELP
   * --version
   * -v --verbose  (which may be repeated for more verbosity)
   * -o --output  (which requires an option argument)
   *
   * the program will have been terminated if unrecognised options are specified,
   * or if an option without GOPT_REPEAT was repeated, or if an option without 
   * GOPT_ARG had an option argument, or if one with GOPT_ARG had no argument.
   */

  if( gopt( options, 'h' ) ){
    /*
     * if any of the help options was specified
     */
    fprintf( stdout, "help text\n" );
    exit( EXIT_SUCCESS );
  }

  if( gopt( options, 'z' ) ){
    /*
     * if --version was specified
     * NB: 'z' is just a key within the source code
     * if you specify -z, it will be treated as an unknown option
     */
    fprintf( stdout, "version number\n" );
    exit( EXIT_SUCCESS );
  }

  if( gopt_arg( options, 'o', & filename ) && strcmp( filename, "-" ) ){
    /*
     * if -o or --output was specified, and its argument was not "-"
     */
    out_stream= fopen( filename, "wb" );
    if( ! out_stream ){
      fprintf( stderr, "%s: %s: could not open file for output\n", argv[0], filename );
      exit( EXIT_FAILURE);
    }
  }
  else
    out_stream= stdout;

  verbosity= gopt( options, 'v' );
  /*
   * return value is the number of times that the option was specified
   */

  if( verbosity > 1 )
    fprintf( stderr, "being really verbose\n" );

  else if( verbosity )
    fprintf( stderr, "being verbose\n" );

  gopt_free( options );
  /*
   * release memory used
   * you can no longer call gopt() etc.
   */

  for( i= 0; i < argc; ++i )
    fprintf( out_stream, "%s\n", argv[i] );
    /*
     * all the options have been removed
     * leaving only the program name and operands
     */

  exit( EXIT_SUCCESS );
}
/* To accept an option that can be repeated and takes an argument, use
 * GOPT_REPEAT | GOPT_ARG.  If you do this then you will want to use
 * gopt_arg_i() or gopt_args().  See gopt.h for details.
 */
