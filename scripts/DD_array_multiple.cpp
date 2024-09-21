// This is a script to calculate a 2D undirected graph containing a full tree, having the posittion of the data points as input.
// output is a two dimensional n*n array, where each row consist of the distance of a particle with every other particle.
// since the distance between the particle and itself is zero, the output array is traceless.
// since we don't want to have duplication (between paricle i and j and j and i), the matrix is upper triangular.
// this script supposed to work with multiple files
// run the script with g++ -std=c++0x DD_array_multiple.cpp

# include <iostream>
# include <fstream>
# include <string>
# include <sstream>
# include <cmath>
# include <cstdlib>
# include <string>

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
  double **DD;
  
  string num, kk;
  string input_filename_base, input_filename_path, input_filename;
  string output_filename_base, output_filename_path, output_filename;
  
  input_filename_path = "/home/nassim/analysis/clustering/position_less40/nope_pi_sn";
  input_filename_base = input_filename_path + "/position_less40_";
  output_filename_path = "/home/nassim/analysis/clustering/dd_less40/nope_pi_sn";
  output_filename_base = output_filename_path + "/dd_less40_";
  for (k = 9; k < 100; k++)
  {
    l = k*10;
    kk = to_string(l);
    if (l < 10)
    {
      num = "00" + kk;
    }
    else if (l < 100)
    {
      num = "0" + kk;
    }
    else if (l < 1000)
    {
      num = kk;
    }
    
    input_filename = input_filename_base + num + ".txt";
    output_filename = output_filename_base + num + ".txt";
    ifstream input (input_filename);
    
      if (input.is_open())
      {
	cout << "reading the input file -- > snapshot:" << kk << "\n";
	
	input >> Nstar;
	input >> maxx;
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
      
      DD = new double*[Nstar];
      for (int i = 0; i < Nstar; ++i)
      {
	DD[i] = new double[Nstar];
      }
      for (i = 0; i < Nstar; i++)
      {
	for (j = 0; j < Nstar; j++)
	{
	  DD[i][j] = 0;
	}
      }
      
      for (i = 0; i < Nstar; i++)
      {
	for (j = i; j < Nstar; j++)
	{
	  DD[i][j] = sqrt((Dx[i] - Dx[j])*(Dx[i] - Dx[j]) + (Dy[i] - Dy[j])*(Dy[i] - Dy[j]));
//	  cout << "DD[" << i << "][" << j << "] = " << DD[i][j] << "\t"; 
	}
      }  
      cout << endl;
      cout << "writing the output into a text file" << endl;
      
	ofstream output (output_filename);
      if (output.is_open())
      {
	cout << "successfully opened the output file" << output_filename << endl;
	for (i = 0; i < Nstar; i++)
	{
	  for (j = 0; j < Nstar; j++)
	  {
	    output << DD[i][j] << " "; 
	  }
	  output << endl;
	}
      output.close();
      }
      else
      {
	cout << "could not open the output file";
      }
      delete[] Dx;
      delete[] Dy;
      delete[] DD;
      
  }
  return 0;
}
      
      
      
      
    
