package fbi.genome.sequencing.rnaseq.simulation;

import java.util.Random;                                    
import java.lang.Math;                                      
import java.util.Arrays;
import java.io.*;     

/* usage 
 * 
 * java Rna_fragmentator pippo.txt 100 0.0 3 > frgs.txt
 *
 *  in order: 
 *
 *  file containing the molecule (pippo.txt)
 *  number of uniform cuts (100)
 *  d0 (0.0)
 *  delta (3)
 *
 *rna_fragmentator
 *----------------   
 *  
 * this program fragments a linear sequence assuming that some 
 * geometrical quantity related to the linear sequence is cut uniformly,
 * and computing the resulting linear fragments (which follow 
 * a Weibull distribution with parameters d0,eta,delta ).
 *
 * see for example
 *
 * A probability concept about size distribution of 
 * sonicated lipid vesicles
 * 
 * (Tenchov,Yanev,Tihova,Koynova) 1985
 *
 *
 *  
 *
 * the observed quantity is called d. The unobservable quantity
 * (which is uniformly cut) is called x.
 * The relation between x and d is given by x = c*(d-d0)^delta
 *
 * eta d0 and delta are parameters of the final Weibull distribution.
 * 
 * the scale parameter of the distribution, eta, and the number of breakpoints
 * are related to c via
 *
 * c = 1/(n*eta^delta)
 *
 * c is determined in such a way that the final fragments normalize to
 * the length of the molecule
 *
 * c=(l-(n+1)d0)/(sum of x_i^(1/delta))
 *
 ************************************************************************************
 * Copyright (c) 2011 Emanuele Raineri                                              
 * distributed under  MIT license                                                               
 *                                                                                  
 * http://en.wikipedia.org/wiki/MIT_License
 *
 *
 ************************************************************************************ 
 */

public class Rna_fragmentator2 {

	  public static void main (String[] args) {
	    System.err.printf("RNA fragmentator\n(c)2011 Emanuele Raineri\n");
	    System.err.printf("reading from:%s\n",args[0]);
	    String line="";
	    try {
	      BufferedReader br = new BufferedReader(new FileReader(args[0]));
	      line = br.readLine();
	      /* cuts per unit length on x, not directly on the molecule */
	      int n=Integer.valueOf(args[1]).intValue();
	      /* parameters of the final Weibull distribution : eta, d0, delta */
	      float d0    = Float.valueOf(args[2]).floatValue();
	      float delta = Float.valueOf(args[3]).floatValue(); 
	      //
	      String [] frags = cut_molecule(line,n,d0,delta);
	      print_string_vector(frags);
	      //System.err.printf("c=%g\n",1/(float)(n*Math.pow(eta,delta)));
	      br.close(); 
	    } catch (IOException e){System.err.println("Error: " + e);}
	  }

	  static void print_float_vector(float[] v){
	    for (int i=0;i<v.length;i++){
	      System.err.printf("%g\n",v[i]);
	    }
	  }

	  static void print_string_vector(String[] v){
	    for (int i=0;i<v.length;i++){
	      if (v[i]!=null && v[i].length()>0) System.out.printf("%s\n",v[i]);
	    }
	  }
	    /**
	      cuts uniformly some geometrical quantity related to the linear
	      dimension of the molecule
	      @param n int number of fragments
	      @return float[]
	    */
	  static float[] uniform_cuts(int n){
	    Random generator = new Random();
	    /* =x= stores the coordinates of the breakpoints, 
	    *  =fragments= stores the lengths of the fragments*/ 
	    float [] fragments = new float [n+1], x = new float [n];
	    /* compute =n= random numbers (uniform distribution) */
	    for(int i=0;i<n;i++){
	      x[i]=generator.nextFloat();  
	    } 
	    /* use them to generate fragments */
	    Arrays.sort(x);
	    /* first fragment */
	    fragments[0]=x[0];
	    for (int i=0;i<n-1;i++){
	      fragments[i+1]= x[i+1]-x[i];
	    }
	    /* last fragment */
	    fragments[n]=1f-x[n-1]; 
	    return fragments;
	  }
	  /**
	 * computes  c so that the fragments of the molecule correctly 
	 * normalize to the molecule length
	 * @param l molecule length
	 * @param n number of breakpoints
	 * @param d0 param of the final weibull
	 * @param x vector of random numbers, uniformly distributed between 0 and 1
	 * @param delta param of the final weibull
	 */
	  static float compute_c (int l, int n, float d0,float [] x, float delta){
	    float sum = 0;
	    for(int i=0;i<x.length;i++){
	      sum+=Math.pow(x[i],1/delta);
	    }  
	    float c = (l - (n+1)*d0)/sum;
	    c=(float)Math.pow(c,-delta);
	    return c;
	  }
	  /**
	  * converts the fragments on the unobservable quantity
	  * into linear fragments of the molecule - which turn out to be Weibull distributed
	  *
	  * @param fragments float array containing the fragments coming from uniform fragmentation
	  * @param d0, delta parameters of the distribution
	  * @param c computed with compute_c
	  * @return float[]
	  *
	  * */
	  static float [] weibull_from_uniform(float [] fragments,  float c, float d0, float delta){
	    int n = fragments.length;
	    float [] d = new float [n];
	    for (int i=0;i<n-1;i++){
	      d[i]=  d0 + (float)Math.pow(fragments[i]/c,(1/delta));
	    }
	    return d; 
	  } 
	  /**
	 * performs the cutting on a string representing an RNA molecule
	 *
	 * @param molecule String the string representing the molecule
	 * @param n int number of cuts to be performed
	 * @param c float constant in the relation x = c*(d-d0)^delta
	 * @param d0 float see explanation for c
	 * @param delta float see explanation for c
	 * @return String[] 
	 */
	  static String[] cut_molecule (String molecule, int n, float d0, float delta){
	    int l = molecule.length();
	    System.err.printf("molecule length:%d\n",l);
	    float[] uni_cuts = uniform_cuts (n);
	    float c = compute_c( l, n, d0,uni_cuts,  delta);
	    float eta = (float)Math.pow((1/(c*n)),(1/delta));
	    System.err.printf("c=%g,eta=%g\n",c,eta);
	    float[] weibull_cuts=  weibull_from_uniform (uni_cuts,c,d0,delta);
	    print_float_vector(weibull_cuts);
	    int current_pos=0;
	    String[] wfrags=new String[n+1];
	    for(int i=0;i<n;i++){ /* n is the number of breakpoints */
	      int to = current_pos+(int)Math.round(weibull_cuts[i]) ;
	      System.err.printf("from %d to %d\n",current_pos,to);
	      // in .substring(startIndex,endIndex), 
	      // startIndex is inclusive, endIndex is exclusive
	      wfrags[i]=molecule.substring(current_pos,to);
	      current_pos=to ;
	      if (current_pos>=l) { 
	        System.err.printf("something is wrong, currentpos bigger than length of molecule \n");
	        System.exit(1);
	      }   
	    }
	    /*last fragment*/
	    System.err.printf("from %d to %d\n",current_pos,l-1);
	    wfrags[n]= molecule.substring(current_pos);
	    return wfrags;
	  }
}
