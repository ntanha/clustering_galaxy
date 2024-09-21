// this is a script to calculate 2 dimensional two-point function with logarithmic bins

# include <iostream>
# include <fstream>
# include <string>
# include <sstream>
# include <cmath>
# include <cstdlib>

using namespace std;


int main()
{
  // reading into an input ascii file, assigning the data to variables.
  // number of stars, size of the box, position of the stars
  
  int i, j, k, l;
  int Nstar;
  double maxx, minx, maxy, miny;  
  
  double *Dx;
  double *Dy;
  
  ifstream input ("./positions/noPE_PI_SN/position_990_17_900to1000_2D.txt");
  
  if (input.is_open())							//check if the file stream was successful in opening a file
  {
    input >> Nstar;							//number of star particles
    input >> maxx;							//boundaries of the box 
    input >> minx;
    input >> maxy;
    input >> miny;
    
    Dx = new double[Nstar];
    Dy = new double[Nstar];
    
    for (i = 0; i < Nstar; i++)
    {
      input >> Dx[i]; 
    }
    for (i = 0; i < Nstar; i++)
    {
      input >> Dy[i];
    }
    
    input.close();
  }
  else
  {
    cout << "could not open the input file";				//give an error message if opening the file was not successful
  }
  cout << "Nstar = " << Nstar <<  ",\t" << minx << " < x < " << maxx << ",\t" << miny << " < y < " << maxy << endl;
   
   
   /* making DD array and fill it with every star-star seperation and binning them
   in binning we want to make sure that the smallest bin size is 2 pc. (From softening length in the simulation).
   to make the code cleaner, I make a map from my real space to log space, which makes my bins linear. 
   I am throwing everything under 2 pc in a separate bin (bin zero) but do not keep it in my calculation. Because everything below that is not trust worthy/numerical result.
   */
   
  double maxsep = sqrt((maxx-minx)*(maxx-minx) + (maxy-miny)*(maxy-miny));
  double minsep = 0;

  double xstart = 2;							//starting point for binning my space
  double xend = maxsep;							//ending point for binning my space
  double logxstart = log10(xstart);				
  double logxend = log10(xend);
  double Lmin = 2;							//minimum length of the smallest bin
  int bins;								//number of bins 
  bins = (int)((log10((xend)/(xstart)))/(log10(1+(Lmin)/(xstart)))) + 1;//the number of bins necessary to make the smallest bin to be 2 pc + the first bin(deltax<2pc)
  Lmin = xstart*(pow((xend/xstart),(1.0/(bins-1))) - 1);		//updates the smallest bin length as we round the number of bins
  double logbinsize = (logxend - logxstart)/(bins-1);			//bin size in the logarithmic space
  
  double *bin;
  bin = new double [bins];
  for (i = 0; i < bins; i++) bin[i] = 0.0;
  for (i = 0; i < bins; i++)
  {
    bin[i] = pow(10, (logxstart + i*logbinsize));
    cout << "bin[" << i << "] = " << bin[i] << "\t";
  }
  
  cout << minsep << "< seperation <" << maxsep << "\n" << bins << " bins, with lengths of " << Lmin << "in logarithmic scale" << endl;
  for(i = 0; i < bins; i++) cout << pow(10, (logxstart + i*logbinsize)) << ", ";
  cout << "\n";

  // the number of seperations in each bin
  double *DD_dis;
  DD_dis = new double [bins];
  for (i = 0; i < bins; i++) DD_dis[i] = 0.0;
  
  double DD;								// distance between each two pair of star particles
  double logDD;
  int DD_index;
  double TotDD = 0;

  
  cout << "calculating the DD array" << endl;
  
  for (i = 0; i < Nstar; i++)
  {
    for (j = 0; j < Nstar; j++)
    {
      DD = sqrt((Dx[i] - Dx[j])*(Dx[i] - Dx[j]) + (Dy[i] - Dy[j])*(Dy[i] - Dy[j]));
      if (DD > xstart)
      {
	logDD = log10(DD);
	DD_index = ((int)((logDD - logxstart)/logbinsize)) + 1;
      }
      else
      {
	logDD = 0;
	DD_index = 0;
      }
      DD_dis[DD_index] = DD_dis[DD_index] + 1;
    }
  }
  
  //normalizing the DD distribution 
  cout << "binning and normalizing the DD array" << endl;
  
  for(i = 0; i < bins; i++)
  {
   
    TotDD = TotDD + DD_dis[i];
  }
  
  for(i = 0; i < bins; i++) DD_dis[i] = DD_dis[i] / TotDD;
  
  for(i = 0; i < bins; i++) cout << "DD[" << i << "] = "<< DD_dis[i] << "\t";
  
  cout << "\n";
  
  //making a random field of postitions of Nstar particles in a box of the size of max and min of star distribution
  
  double DR;							// Distance between each two star particles and the random field 
  double RR;							// Distance between each two pair of random field
  double logDR;
  double logRR;
  
  //postitons of the random field
  double *Rx;
  Rx = new double[Nstar];
  
  double *Ry;
  Ry = new double[Nstar];
  
  //the number of seperations in each bin
  double *RR_dis;
  RR_dis = new double [bins];
  for (i = 0; i < bins; i++) RR_dis[i] = 0.0;
  
  double *DR_dis;
  DR_dis = new double [bins];
  for (i = 0; i < bins; i++) DR_dis[i] = 0.0;
  
  double TotRR = 0;
  double TotDR = 0;
  int RR_index;
  int DR_index;  
  
  int num = 1000;						//number of random fields to average on
  
  cout << "calculating the RR and DR array" << endl;
  
  // calculating the position of the random star particles
  for (i = 0; i < num; i++)
  {
    for (j = 0; j < Nstar; j++)
    {
      Rx[j] = ((double)(rand() % RAND_MAX) / RAND_MAX) * (maxx-minx) + minx;
      Ry[j] = ((double)(rand() % RAND_MAX) / RAND_MAX) * (maxy-miny) + miny;
    }
    
    // calculating the distances between every two random star particles and between the random field and star particles
    for (j = 0; j < Nstar; j++)
    {
      for (k = 0; k < Nstar; k++)
      {
	RR = sqrt((Rx[j] - Rx[k])*(Rx[j] - Rx[k]) + (Ry[j] - Ry[k])*(Ry[j] - Ry[k]));
	
	DR = sqrt((Dx[j] - Rx[k])*(Dx[j] - Rx[k]) + (Dy[j] - Ry[k])*(Dy[j] - Ry[k]));
	
	if (RR > xstart)
	  {
	    logRR = log10(RR);
	    RR_index = ((int)((logRR - logxstart)/logbinsize)) + 1;
	  }
	else
	  {
	    logRR = 0;
	    RR_index = 0;
	  }
      	if (DR > xstart)
	  {
	    logDR = log10(DR);
	    DR_index = ((int)((logDR - logxstart)/logbinsize)) + 1;
	  }
	else
	  {
	    logDR = 0;
	    DR_index = 0;
	  }
	
	//binning the RR
	RR_dis[RR_index] = RR_dis[RR_index] + 1;
	DR_dis[DR_index] = DR_dis[DR_index] + 1;
	
      }
    }
  }
  
  //averaging over the different random fields and normalizing the RR field.
  
  cout << "averaging and normalizing RR and DR arrays" << endl;
  
  for (i = 0; i < bins; i++)
  { 
    TotRR = TotRR + RR_dis[i];
    TotDR = TotDR + DR_dis[i];
  }
  
  for (i = 0; i < bins; i++)
  {
    RR_dis[i] = RR_dis[i] / TotRR;
    DR_dis[i] = DR_dis[i] / TotDR;
  }
  
  for(i = 0; i < bins; i++) cout << "RR[" << i << "] = "<< RR_dis[i] << "\t";
  cout << "\n";
  for(i = 0; i < bins; i++) cout << "DR[" << i << "] = "<< DR_dis[i] << "\t";
  cout << "\n";
  
  //Calculating the two-point function
  
  cout << "calculating the twopoint function" << endl;
  
  double *TwoPoint;
  TwoPoint = new double[bins];
  for (i = 0; i < bins; i++) TwoPoint[i] = 0.0;
  
  for (i = 0; i < bins; i++)
  {
    if (RR_dis[i] != 0)
      {
	TwoPoint[i] = (DD_dis[i] - 2*DR_dis[i] + RR_dis[i])/RR_dis[i];
      }
  }
  for(i = 0; i < bins; i++) cout << "TwoPoint[" << i << "] = "<< TwoPoint[i] << "\t";
  cout << "\n";
  
  //writing the output into a text file
  
  cout << "writing the output into a text file" << endl;
  
    ofstream output ("./2point_data/noPE_PI_SN/2point_990_log_2D_900to1000.txt");
  if (output.is_open())
  {
    for (i = 0; i < bins; i++)
    {
      output << bin[i] << " " << DD_dis[i] << " " << DR_dis[i] << " " << RR_dis[i] << " " << TwoPoint[i] << endl;
    }
  output.close();
  }
  else
  {
    cout << "could not open the output file";
  }
  
  return 0;
}
