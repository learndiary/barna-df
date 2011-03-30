package fbi.genome.sequencing.rnaseq.simulation;

import java.util.Random;                                    
import java.lang.Math;                                      
import java.util.Arrays;
import java.io.*;     

/* usage 
 * 
 *  java Rna_fragmentator pippo.txt 10000 1.0 0.0 3
 *
 *  in order: 
 *
 *  file containing the molecule (pippo.txt)
 *  number of uniform cuts (10000)
 *  c (1.0)
 *  d0 (0.0)
 *  delta (3)
 *
 *
 * */


/**
 * this program fragments a linear sequence assuming that some 
 * geometrical quantity related to the linear sequence is cut uniformly,
 * and computing the resulting linear fragments (which follow 
 * a Weibull distribution with parameters d0,eta,delta ,defined as below).
 *
 * see for example
 *
 * A probability concept about size distribution of 
 * sonicated lipid vesicles
 * 
 * (Tenchov,Yanev,Tihova,Koynova) 1985
 *
 *
 * the MIT license 
 *
 * Copyright (c) 2011 Emanuele Raineri, MIT license 
 *
 *
 * 
 *
 * the observed quantity is called d. The unobservable quantity
 * (which is uniformly cut) is called x.
 * The relation between x and d is given by x = c*(d-d0)^delta
 *
 * d0 and delta are parameters of the final Weibull distribution.
 * The scale parameter of the Weibull, eta is given by 
 * 
 * eta=L*(1/(c*n))^(1/delta)
 *
 * where n is the number of cuts, L is the length of the molecule
 * 
 */

public class Rna_fragmentator {

  public static void main (String[] args) {
    System.err.printf("RNA fragmentator\nby Emanuele Raineri\n");
    System.err.printf("reading from:%s\n",args[0]);
    String line="";
    try {
      BufferedReader br = new BufferedReader(new FileReader(args[0]));
      line = br.readLine();
      int numberOfCuts=Integer.valueOf(args[1]).intValue();
      float c = Float.valueOf(args[2]).floatValue();
      float d0= Float.valueOf(args[3]).floatValue();
      float delta = Float.valueOf(args[4]).floatValue(); 
      String [] frags = cut_molecule(line,numberOfCuts,c,d0,delta);
      print_string_vector(frags);
      System.err.printf("eta=%g\n",line.length()*Math.pow((1/(c*numberOfCuts)),(1/delta))); 
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
      @param minvol float smallest possible fragment
      @param vol total volume
      @param n int number of fragments
      @return float[]
    */
  static float[] uniform_cuts(float minvol, float vol, int n){
    Random generator = new Random();
    /* =x= stores the coordinates of the breakpoints, 
    *  =fragments= stores the lengths of the fragments*/ 
    float [] fragments = new float [n+1], x = new float [n];
    /* compute =n= random numbers (uniform distribution) */
    float deltavol=(vol-minvol);
    for(int i=0;i<n;i++){
      x[i]=generator.nextFloat()*deltavol+minvol;  
    } 
    /* use them to generate fragments */
    Arrays.sort(x);
    /* first fragment */
    fragments[0]=x[0];
    for (int i=0;i<n-1;i++){
      fragments[i+1]= x[i+1]-x[i];
    }
    /* last fragment */
    fragments[n]=vol-x[n-1]; 
    return fragments;
  }

  /**
  * converts the fragments on the unobservable quantity
  * into linear fragments of the molecule - which turn out to be Weibull distributed
  *
  * @param fragments float array containing the fragments coming from uniform fragmentation
  * @param l int length of the molecule
  * @return float[]
  *
  * */
  static float [] weibull_from_uniform(float [] fragments, float c, float d0, float delta){
    int n = fragments.length;
    float [] d = new float [n];
    double sum= 0d;
    for (int i=0;i<n-1;i++){
      d[i]=  d0 + (float)Math.pow(fragments[i]/c,(1/delta));
      sum+= d[i];
    }
    System.err.println("sum of linear lengths: "+ sum);
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
  static String[] cut_molecule (String molecule, int n, float c,float d0, float delta){
    int l = molecule.length();
    System.err.printf("molecule length:%d\n",l);
    float vol =  c*(float)Math.pow(l,delta);
    float minvol = 0; // c*(float)Math.pow(d0,delta);
    float[] uni_cuts = uniform_cuts (minvol,vol,n);
    float[] weibull_cuts=  weibull_from_uniform (uni_cuts,c,d0,delta);
    //print_float_vector(weibull_cuts);
    int current_pos=0;
    String[] wfrags=new String[n];
    for(int i=0;i<n;i++){ /* n is the number of breakpoints */
      int to = current_pos+(int)weibull_cuts[i];
      if (to>=l) {
          /*System.out.printf("overflow!\n");*/
          to=l-1;
      } 
      //System.out.printf("from %d to %d\n",current_pos,to);
      wfrags[i]=molecule.substring(current_pos,to);
      current_pos=current_pos+(int)weibull_cuts[i];
      if (current_pos>=l) break;  
    }
    /*last fragment*/
    //wfrags[n]= molecule.substring(current_pos,l+1);
    return wfrags;
  }
}
