#include "R.h"
#include "Rmath.h"

#include "Const.h"
#include "Vector.h"
#include "Matrix.h"


int debug = 0;

void gelmanRubin( CMatrix * x,
                  double * gelman_rubin_sqrt_R,
                  double * gelman_rubin_sqrt_upper_quartile )
{
  int i;
  double overall_mean, between_var, within_var, var_within_var, var_between_var, cov_within_between;
  double hat_V, hat_var_V, prop_w, prop_b, prop_wb;
  double degrees_freedom_within_var, degrees_freedom_hat_V;
  double minimumDegreesFreedom = 3;

  //x->Col() is the number of chains
  //x->Row() is the number of samples per chain

  CVector mean( x->Col() );
  CVector var( x->Col() );

  overall_mean = 0;
  for ( i = 0; i < x->Col(); i++ )
  {
    mean.Val(i) = x->getColumn(i).Mean();
    var.Val(i) = x->getColumn(i).Var();
    overall_mean += mean.Val(i);
  }
  overall_mean /= (double) x->Col();

  between_var = 0;
  within_var = 0;
  for ( i = 0; i < x->Col(); i++ )
  {
    between_var += ( mean.Val(i) - overall_mean ) * ( mean.Val(i) - overall_mean );
    within_var += var.Val(i);
  }
  between_var *= ( ( (double) x->Row() ) / ( (double) x->Col() - 1.0 ) );
  within_var /= ( (double) x->Col() );

  //estimate var V (i.e. get hat(var(V)) )
  var_within_var = var.Var() / ( (double) var.Len() );
  
  var_between_var = between_var * between_var * ( 2.0 / ( (double) x->Col() - 1.0 ) );

#ifdef FIX1
  CVector tmpvec = mean.square();
  cov_within_between = ( var.Cov( tmpvec ) - 2 * overall_mean * var.Cov( mean ) );
#else
  cov_within_between = ( var.Cov( mean.square() ) - 2 * overall_mean * var.Cov( mean ) );
#endif
//  cov_within_between = ( var.Cov( mean.square() ) - 2 * overall_mean * var.Cov( mean ) );
  cov_within_between *= ( ( (double) x->Row() ) / ( (double) x->Col() ) );

  if ( debug )
  {
#ifdef FIX1
    CVector tmpvec = mean.square();
    printf( "var.Cov( mean.square() ) =  %f,  var.Cov( mean ) = %f,  overall_mean = %f\n",
      var.Cov( tmpvec ),
#else
    printf( "var.Cov( mean.square() ) =  %f,  var.Cov( mean ) = %f,  overall_mean = %f\n",
      var.Cov( mean.square() ),
#endif
//	    var.Cov( mean.square() ),
      var.Cov( mean ), overall_mean );
  }

  prop_w = ( (double) x->Row() - 1 ) / ( (double) x->Row() );
  prop_b = ( (double) x->Col() + 1 ) / ( (double) x->Col() * x->Row() );
  prop_wb = ( ( (double) x->Col() - 1 ) * ( (double) x->Row() - 1 ) )
          / ( (double) x->Col() * x->Row() * x->Row() );

  hat_var_V = var_within_var * prop_w * prop_w 
            + var_between_var * prop_b * prop_b 
            + 2 * prop_wb * cov_within_between;

  if ( debug )
  {
    printf( "var_within var = %f,  var_between var = %f,  cov_within_between = %f\n",
	    var_within_var, var_between_var, cov_within_between );
  }

  //estimate degrees of freedom for within var
  if ( within_var * within_var >= 0.5 * DELTA * var_within_var )
  {
    degrees_freedom_within_var = 2.0 * within_var * within_var / var_within_var;
  }
  else 
  {
    printf( "gelmanRubin: degrees of freedom of within variance W estimate is too small. Set to 1.0.\n" );
    degrees_freedom_within_var = 1.0;   
  }

  //estimate degrees of freedom for hat(V)/V
  hat_V = ( within_var * prop_w )  + ( between_var * prop_b );
  if ( hat_V * hat_V >= 0.5 * DELTA * hat_var_V )
  {
    degrees_freedom_hat_V = 2.0 * hat_V * hat_V / hat_var_V;
  }
  else
  {
    printf( "gelmanRubin: degrees of freedom of V estimate is too small.\n" );
    degrees_freedom_hat_V = DELTA;
  }

  if ( debug )
  {
    printf( " hat_V = %f, hat_var_V = %f, within_var = %f, between_var = %f, dfhatV = %f\n",
	    hat_V, hat_var_V, within_var, between_var, degrees_freedom_hat_V );
  }

  //make sure degrees of freedom is larger than 2 (use correction)
  degrees_freedom_hat_V += minimumDegreesFreedom;

  (*gelman_rubin_sqrt_R) = ( hat_V / within_var ) 
                         * ( degrees_freedom_hat_V / ( degrees_freedom_hat_V - 2.0 ) );
  (*gelman_rubin_sqrt_R) = sqrt( (*gelman_rubin_sqrt_R) );

  (*gelman_rubin_sqrt_upper_quartile) = prop_w 
                 + prop_b * qf( 0.975, (double) x->Col() - 1, degrees_freedom_within_var, 1, 0 );
  (*gelman_rubin_sqrt_upper_quartile) *= ( degrees_freedom_hat_V / ( degrees_freedom_hat_V - 2.0 ) );
  (*gelman_rubin_sqrt_upper_quartile) = sqrt( (*gelman_rubin_sqrt_upper_quartile) );

  if ( debug )
  {
    printf("\nprop_b = %f, prop_w = %f, prop_wb = %f, dfhatV = %f dfwv = %f\n",
	   prop_b,prop_w,prop_wb,degrees_freedom_hat_V,degrees_freedom_within_var );
    printf("gelman_rubin_sqrt_R = %f,  gelman_rubin_sqrt_upper_quartile = %f\n",
	   (*gelman_rubin_sqrt_R), (*gelman_rubin_sqrt_upper_quartile) );
  }

}//end gelmanRubin




void gelmanRubinDiagnostics( double * data,
                             long * number_chains,
                             long * number_samples,
                             long * sample_period,
                             double * gelman_rubin_sqrt_R,
                             double * gelman_rubin_sqrt_upper_quartile )
{
  bool finished;
  int count, n_samples, starting_index, ending_index;
  double sqrt_R, upper_quartile;
  CMatrix * sub_x;

  CMatrix x( data, (int) (*number_samples), (int) (*number_chains), false );

  if ( debug > 2 )
  {
    printf( "gelmanRubinDiagnostics: data matrix is: \n" );
    x.Print();
  }

  starting_index = 0;
  ending_index = (int) (*sample_period);

  count = 0;
  finished = false;
  while ( !finished )
  {
    //get first rows of data
    n_samples = ending_index - starting_index;
    sub_x = new CMatrix( x.getRows(0, ending_index - 1 ) );

    //compute the Gelman-Rubin statistics
    gelmanRubin( sub_x, &sqrt_R, &upper_quartile );
    
    //store the results
    gelman_rubin_sqrt_R[ count ] = sqrt_R;
    gelman_rubin_sqrt_upper_quartile[ count ] = upper_quartile;

    if ( debug )
    {
      printf("gelman_rubin_sqrt_R[%d] = %f,  gelman_rubin_sqrt_upper_quartile = %f\n",
	     count, gelman_rubin_sqrt_R[count], gelman_rubin_sqrt_upper_quartile[count] );
    }

    count++;

    delete sub_x;

    //check if another statistics can be computed
    ending_index += (int) (*sample_period);
    if ( ending_index > (int) (*number_samples) )
    {
      finished = true;

      if ( debug )
      {
        printf( "count = %d\n", count);
      }
    }
  }//end loop

}//end

