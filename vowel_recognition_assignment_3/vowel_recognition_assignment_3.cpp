// vowel_recognition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define P 12
#define N 320
#define framecount 5
char vowel[]={'a','e','i','o','u','\0'};
	  
long double PI =  3.1428571428571;// Global value of PI.
double weights[] = {1.0,3.0,7.0,13.0,19.0,22.0,25.0,33.0,42.0,50.0,56.0,61.0};
int StartIndex =0;
//here we will set start index of each recording
void findStart(FILE* filepointer)
{   long double max=0.0;
	StartIndex=0;
	char input[20];
	long double samplevalue=0.0;
	if(!filepointer)
		{printf("start index cannnot be found\n");
	       return ;
	      }
	
		   int i=10,count=0;
		  
		/*
		we are using i to skip 10 strings as they are discription and not data 
		then we are using framesize to count upto 320 samples in each frame 
		we are using count to store count of frames when max energy moved from its previous frame to new frame 
		*/
		  while(!feof(filepointer))
	     {  
		   fscanf(filepointer,"%s",input);
		   samplevalue=atof(input);
		   if(i>0)
		   {
			   i--;
			  
		   }
		   else
		   {
			      if(max < samplevalue)
				      { 

						  max=samplevalue;
				        
						  StartIndex+=count;
				          count=0;     
				     }
				   else
				     count++;
				  
			   }
            }//end while
	
	// after the function returns we would have startIndex of each recording 
}
//here we will do dc shift on text file and return vlaue of dc shift
long double dcShift(FILE * filepointer)
{   

	long double sampleValue, sampleSum=0.0,sampleCount=0;
   char input[20];
	if(!filepointer)
		{printf("DC shift cannot be done\n");
	       return 0.0;}
	
	else {
		   int i=10;
		//skip first five lines as they are not data
	     while(!feof(filepointer))
	     {  
		   fscanf(filepointer,"%s",input);
		   if(i>0)
		   {
			   i--;
			  
		   }
		   else
		   {
	      //here we will start reading the samples and then find sum of all values and no of samples  
		    
		   sampleCount++;
		   sampleValue=atof(input);
		   sampleSum+=sampleValue;
		   }
	     }
	  }
	fclose(filepointer);	
	long double dcOffset=sampleSum/sampleCount;
		return dcOffset;
}
//  normalization  on text file and return value of normalization factor

 long double normalize(FILE * filepointer)
{  
	
	long double sampleValue,sampleCount=0,max=0.0,min=0.0;
	char input[20];

	if(!filepointer)
		printf("normaization cannot be done\n");
	    
	else {
		   int i=10;
		//skip first five lines as they are not data
	   while(!feof(filepointer))
	   {  
		   fscanf(filepointer,"%s",input);
		   if(i>0)
		   {
			   i--;
			  
		   }
		   else{
	      //here we will start reading the samples and then find sum of all values and also find min and max values   
		    
		   sampleValue=atof(input);
		   if(max< sampleValue)
			  max=sampleValue;
		   if(min> sampleValue)
			  min=sampleValue;
		   }
	   } 
	  min=abs(min);
	long double avg=(long double)(max+min)/2;
	 long double normalizeFactor=(long double)((max-5000)/max);
	//normalizeFactor=(int) (normalizeFactor*100+0.5);
	//normalizeFactor=(float)normalizeFactor/100;//rounding to one decimal place
	fclose(filepointer);
	return normalizeFactor;

         
     }

 }


 //here i will pass file name and array so that it can store values in array
 void FillArray(FILE * filepointer,long double array[])
 {   
	 //empty the array
	  for(int i=0;i<N;i++)
		 array[i]=0.0;
     
	  char sample[30];
	  long double sampleValue=0.0;
	  int i=0;//read only 320 samples
      while(!feof(filepointer))
			{  //skip tp the start Indexs previous 2
				
	            fscanf(filepointer,"%s",sample);
				 if(i==N)
			         return;
		        	else
						{
							sampleValue=atof(sample);//string to float
						    //now i will start putting 320 values in a array 
				              array[i++]=sampleValue;
					     }     
			}
			
   
	 }

void calculate_Cis(long double arr[],long double getCi[])
{

	//first we will calculate Ris 
    long double Ri[13];
	long double Ai[12];
	long double Ci[13];
	for(int i=0;i<=P;i++)
	{
		long double sum=0.0;
		for(int j=0;j<N-i;j++)
		{
			 sum=sum+arr[j]*arr[j+i];
				
		  }
	     	Ri[i]=sum;
		   sum=0;
	}


//	now we will calculate Ais using lavinson durbins algorithm
	long double E[13],K[13],A[13][13];
	
	//initialize base values
	for(int i=0;i<=P;i++)
	{
		E[0]=0;
		K[i]=0;
		
	   for(int j=0;j<=P;j++)
		   A[i][j]=0;
	}
	
	E[0]=Ri[0];

	//repeat for i =0 to p
	for(int i=1;i<=P;i++)
	{
	   for(int j=1;j<i;j++)
		   K[i]+=(A[i-1][j]*Ri[i-j]);
        
	   K[i]=((Ri[i]-K[i])/E[i-1]);
	   A[i][i]=K[i];
	    
	   for(int j=1;j<i;j++)
		   A[i][j]=(A[i-1][j]-(K[i]*A[i-1][i-j]));
      
	   E[i]=(1-K[i]*K[i])*E[i-1];
	}

	//put value in Ai array
	for(int i=1;i<=P;i++)
		Ai[i-1]=A[12][i];
	
	
	
	//now we calculate Cis 
	for(int i=1;i<=P;i++)
		Ci[i]=0;
	Ci[0]=log(Ri[0]*Ri[0]);

	for(int i=1;i<=P;i++)
	{ 
		for(int j=1;j<i;j++)
			Ci[i]+=((double)j/i)*Ci[j]*Ai[i-j-1];

	   Ci[i]+=Ai[i-1];
	   
	}

    
	
	
	//after that i will store these values in a getCi array 
	for(int i=1;i<=P;i++)
		getCi[i-1]=Ci[i];
	
}
//we will get the avg ci of a vowel here 
void getAvgCi(long double ALLCi[10][5][12],long double AvgCi[5][12])
{
	//set avgcgi to 0
	for(int i=0;i<5;i++)
	{
	   for(int j=0;j<P;j++)
	    AvgCi[i][j]=0;
  	}

	//put values in AVgci array
	for(int p=0;p<5;p++)
	{
	  for(int q=0;q<10;q++)
	  {
	    for(int r=0;r<P;r++)
			AvgCi[p][r]+=ALLCi[q][p][r];
	   }
	}
			
	// now i will divide each value by 10 to get avg
			
   for(int i=0;i<5;i++)
	{
	   for(int j=0;j<P;j++)
	    AvgCi[i][j]=AvgCi[i][j]/10;
  	}

}
//fucntion to calculate raised sined window of ci

void RaisedSineWindow(long double AvgCi[5][12])
{
   for(int i=1;i<=12;i++)
   {
   for(int j=0;j<5;j++)
     AvgCi[j][i-1]*=(1+(6*sin((PI*i)/12)));
   }
}

//function to predict the vowel
void findVowel(long double testArray[5][12])
{
	long double vowelPredictArray[5]={0.0,0.0,0.0,0.0,0.0};
	long double refArr[5][12];
	char refFileName[30];
	 long double minValue;
	 int minIndex=0;
	

 
	//open files one by one
	 
	 for(int i=0;i<5;i++)
	{
		
		char input[30];
		long double value=0.0;
		char VowelStr[3]={vowel[i],'\0'};
		sprintf(refFileName,"%s_avgCi.txt",VowelStr);
	    FILE*filepointer=fopen(refFileName,"r");//file opened
	   if(!filepointer)
		{
			printf("file %s \t reFileName \t cannot be opened",refFileName);
	        return;
	   }
	  //fill ref array from refernce file
	    for(int row=0;row<5;row++)
	  {
		  for(int col=0;col<P;col++)
		  {
	        fscanf(filepointer,"%s",input);
		    value=atof(input);
			refArr[row][col]=value;
		  }
	 }
	    //now for all 5 rows calculate tokhura distance
		long double AllDTcep[5]={0.0,0.0,0.0,0.0,0.0};
	for(int row=0;row<5;row++)
	  {
		
		long double DTcep=0.0;
		  for(int col=0;col<P;col++)
		  {
			  DTcep+=weights[col]*((testArray[row][col]-refArr[row][col])*(testArray[row][col]-refArr[row][col]));
            }
		  //put this in a array AllDTcep[5];
		  AllDTcep[row]=DTcep/12;
	 }
	 // now take average of AllDTcep array and store it in vowelPredictArray 	  
   long double sum=0.0;
	for(int itr=0;itr<5;itr++)
		sum+=AllDTcep[itr];
	sum=sum/5;

	vowelPredictArray[i]=sum; //i have enterd avg tokhura distance of 5 rows in a variable sum and saved it in array at index i  
	fclose(filepointer);
	}//repeat for all vowel reference files	 
	
	
	//find min value in array and print which vowel is it
    minValue=vowelPredictArray[0];//set minvalue at index 0 
	for(int i=0;i<5;i++)
	{ 
		if(minValue>vowelPredictArray[i])
		{
		  minIndex=i;
		  minValue=vowelPredictArray[i];
		 }
	}
	//now print vowel at index minIndex
	printf("\t pridiction is \t %c",vowel[minIndex]);
	      
}






int _tmain(int argc, _TCHAR* argv[])
{ 
	 
     long double sampleArr[320];
	  long double ALLCi[10][5][12];
	 long double AvgCi[5][12];
	  char fileName[100];
	  char  sample[30];
      long double testArray[5][12];
	  long double getCi[12];

	  int skipvalues=0;
	  //here we  run this loop for each vowel
	FILE * filepointer;
	for(int i=0;i<5;i++)//for vowel a to u
	{
		char vowelString[3]={vowel[i],'\0'};
		//each vowel has 10 recordings
	
		for(int j=1;j<=10;j++)
		{  

			//here we open 50 files one by one and calculate Ci matrix of each vowel for 10 recordings and store the matrix in 5 separate files 
		     
			sprintf(fileName,"vowels\\214101063_%s_%d.txt",vowelString,j);
			filepointer=fopen(fileName,"r");
			if(!filepointer)
			{	printf("error while opening file\n");
			     break;
			}
			//we first call the function findStart to get the start sample value
			findStart(filepointer);
			
			
			rewind(filepointer);//reset file pointer to start of file
			//here we will call normalize function which will creater new files with normalized values
		    long double dcOffset= dcShift(filepointer);
			//initialize file pointer again
			filepointer=fopen(fileName,"r");
			long double normalizeFactor=normalize(filepointer);
			filepointer=fopen(fileName,"r");// here after normalization function call the file pointer needs to be set again at start of file  
		    		
		  //here i will set the file pointer of this recording according to startIndex
	
	        skipvalues= (StartIndex)+10;//we will move our data till skipvalues no of samples
	      
			//this is to skip first 10 strings 
       
			while(!feof(filepointer))
			{  
			    fscanf(filepointer,"%s",sample);
				if(skipvalues)
				{
				   skipvalues--; 
			      }
				else
				{
					break;
				}
	
			}	
			
			//now my file pointer is at skip values from start
		   //repeat for 5 frames
			for(int k=0;k<5;k++)
		     {  
				  //empty the sample array
				 for(int s=0;s<320;s++)
					 sampleArr[s]=0.0;
				  // fill values in sample array
			     	FillArray(filepointer,sampleArr);
			    
					 //now normalize the sample array 

				for(int i=0;i<N;i++)
				{
				    sampleArr[i]=(sampleArr[i]-dcOffset)*normalizeFactor;
				}
				

				//here i will call calculate Ci function and it will fill the getCi array with required Ci value now i will store this on my AllCi array 
				calculate_Cis(sampleArr,getCi);
			  //do raised sine window here on getCi
				for(int m=1;m<=P;m++)
				{  
					long double s=m;
					getCi[m-1]*=(1+(6*sin((PI*s)/12)));
				}
				//put getCi array value in ALLci array
				
				for(int m=0;m<P;m++)
				ALLCi[j-1][k][m]=getCi[m];
			}// 5 frames per recording done
		 }//one vowel done
	
		
			
         getAvgCi(ALLCi,AvgCi);//calculate avg of 50 rows of ci in 5 rows

		   

		    char output[40];
			sprintf(output,"%s_avgCi.txt",vowelString);
			FILE*outfile=fopen(output,"w");

		 //writing in a text file 
		 
			for(int i=0;i<5;i++)
			{
			   for(int j=0;j<P;j++)
			   {
			      fprintf(outfile,"%lf ",AvgCi[i][j]);
			   }
			   fprintf(outfile,"\n");
			}
			printf("training of %s \t is done\n",vowelString);
			fclose(outfile);
          
}//5 vowels done 
	 
	  
	 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//here comes testing part where i will take one recording at a time and find its Ci values then i will find its distance using tokhura distance 

	//  with respect to all vowels Ci values.
	// then which ever file gives lowest distance is out vowel


	for(int i=0;i<5;i++)//no of vowels
	{  
	   char vowelString[3]={vowel[i],'\0'};
	   for(int j=11;j<=20;j++)
	   
	   { sprintf(fileName,"vowels\\214101063_%s_%d.txt",vowelString,j);
			filepointer=fopen(fileName,"r");
			if(!filepointer)
			{	printf("error while opening file\n");
			     break;
			}
		
			

			findStart(filepointer);// it will set startIndex value
			
			filepointer=fopen(fileName,"r");//reset file pointer to start of file
			//here we will call normalize function which will creater new files with normalized values
		    long double dcOffset= dcShift(filepointer);
			//initialize file pointer again
			filepointer=fopen(fileName,"r");
			long double normalizeFactor=normalize(filepointer);
			filepointer=fopen(fileName,"r");// here after normalization function call the file pointer needs to be set again at start of file  
		    		
		  //here i will set the file pointer of this recording according to startIndex
	
	        skipvalues= StartIndex+10;//we will move our data till skipvalues no of samples
	         
			
			while(!feof(filepointer))
			{  
			    fscanf(filepointer,"%s",sample);
				if(skipvalues)
				{
				   skipvalues--; 
			      }
				else
				{
					break;
				}
	
			}	
			
			//now my file pointer is at skip values from start
		   //repeat for 5 frames
			for(int k=0;k<5;k++)
		     {  
				  //empty the sample array
				 for(int s=0;s<320;s++)
					 sampleArr[s]=0.0;
				  // fill values in sample array
			     	FillArray(filepointer,sampleArr);
			    
					 //now normalize the sample array 

				for(int i=0;i<N;i++)
				{
				    sampleArr[i]=(sampleArr[i]-dcOffset)*normalizeFactor;
				}
				

				//here i will call calculate Ci function and it will fill the getCi array with required Ci value 
				calculate_Cis(sampleArr,getCi);
			   
			//now put this array value in a array of value 5*12
				for(int m=0;m<P;m++)
					testArray[k][m]=getCi[m];
				
			}//one recording done i got the 5 by 12 array of ci values
			
			//apply raised sined window
			 RaisedSineWindow(testArray);
		    
			 // call function which will find tokohura distance
			 printf("\nvowel being tested is\t %s",vowelString);
			 
			 findVowel(testArray);
			
	   }//one vowel all recordings tested
	 
	}//all vowels taken 

	
return 0;
}


