import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;


public class Overlap {

	static class IAcomparator implements Comparator {
		public int compare(Object o1, Object o2) {
			int[] a= (int[]) o1, b= (int[]) o2;
			if (a[0]< b[0]) {
				return -1;
			} else if (a[0]== b[0]) {
				if (a[1]< b[1])
					return -1;
				else if (a[1]> b[1])
					return 1;
				return 0;
			}
			return 1;
		}
	}
	
	static class Startcomparator implements Comparator {
		public int compare(Object o1, Object o2) {
			int[] a= (int[]) o1, b= (int[]) o2;
			if (a[0]< b[0]) {
				return -1;
			} else if (a[0]== b[0]) {				
				return 0;
			}
			return 1;
		}
	}
	public static void overlapPositionToBed(File bed, File pos) {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(pos));
			HashMap<String, List> map= new HashMap<String, List>();
			for (String s = null; (s= buffy.readLine())!= null; ) {
				int p= s.indexOf(':'), q= s.indexOf('-');
				String chr= s.substring(0, p);
				List<int[]> v= map.get(chr);
				if (v== null) {
					v= new Vector<int[]>();
					map.put(chr, v);
				}
				int[] a= new int[3];
				a[0]= Integer.parseInt(s.substring(p+1,q));
				a[1]= Integer.parseInt(s.substring(q+1,s.length()));
				v.add(a);
			}
			buffy.close();
			
			Iterator<List> ii= map.values().iterator();
			Comparator c= new IAcomparator();
			while(ii.hasNext())
				Collections.sort(ii.next(), c);
			
			buffy= new BufferedReader(new FileReader(bed));
			c= new Startcomparator();
			HashMap<String, Integer> map2= new HashMap<String, Integer>();
			for (String s = null; (s= buffy.readLine())!= null; ) {
				String[] ss= s.split("\\s");
				Vector<int[]> v= (Vector<int[]>) map.get(ss[0]);
				if (v== null)
					continue;
				int[] a= new int[3];
				a[0]= Integer.parseInt(ss[1]);
				a[1]= Integer.parseInt(ss[2]);
				int p= Collections.binarySearch(v, a, c);
				if (p< 0)
					p= -(p+1);
				String id= null;
				for (int i = p; i < v.size(); ++i) {
					int[] b= v.get(i);	// probe
					id= null;
					if (a[0]>b[1]|| a[1]<b[0])
						break;
					if (a[0]<= b[1]&& b[0]<=a[1]) {
						String[] start= ss[11].split(",");
						String[] len= ss[10].split(",");
						boolean stop= false;
						for (int j = 0; j < len.length&& !stop; j++) {
							int lenx= Integer.parseInt(len[j]);
							int startx= Integer.parseInt(start[j]);
							for (int k = 0; k < lenx; k++) {
								int x= a[0]+ startx+ k;
								if (x>= b[0]&& x<= b[1]) {
									if (id== null)
										id= ss[0]+":"+b[0]+"-"+b[1]; 
									Integer ix= 0;
									if (map2.containsKey(id))
										ix= map2.get(id);
									map2.put(id, ix+1);
									stop= true;	// hit only once
									break;
								}
								
							}
						}
					}
				}
				
			}
			buffy.close();
			
			Iterator<String> iter= map2.keySet().iterator();
			while(iter.hasNext()) {
				String id= iter.next();
				System.out.println(id+"\t"+map2.get(id));
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		
		args=new String[] { 
//				"/Users/micha/projects/Esterller/H.sapiens/hg19_cpgIslandExt.bed",
//				"/Users/micha/projects/Esterller/H.sapiens/genexpr/all_promoters.pos"};
		
				"/Users/micha/projects/Esterller/M.musculus/mm9_cpgIslandExt.bed",
				"/Users/micha/projects/Esterller/M.musculus/genexpr/all_promotors.pos" 
		};
		
		overlapPositionToBed2(new File(args[0]), new File(args[1]));
	}

	// takes the IDs, not occcurences
	public static void overlapPositionToBed2(File bed, File pos) {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(pos));
			HashMap<String, List> map= new HashMap<String, List>();
			for (String s = null; (s= buffy.readLine())!= null; ) {
				int p= s.indexOf(':'), q= s.lastIndexOf('-');
				String chr= s.substring(0, p);
				List<int[]> v= map.get(chr);
				if (v== null) {
					v= new Vector<int[]>();
					map.put(chr, v);
				}
				int[] a= new int[3];
				a[0]= Integer.parseInt(s.substring(p+1,q));
				a[1]= Integer.parseInt(s.substring(q+1,s.length()));
				v.add(a);
			}
			buffy.close();
			
			Iterator<List> ii= map.values().iterator();
			Comparator c= new IAcomparator();
			while(ii.hasNext())
				Collections.sort(ii.next(), c);
			
			buffy= new BufferedReader(new FileReader(bed));
			c= new Startcomparator();
			HashMap<String, String> map2= new HashMap<String, String>();
			for (String s = null; (s= buffy.readLine())!= null; ) {
				String[] ss= s.split("\\s");
				Vector<int[]> v= (Vector<int[]>) map.get(ss[0]);
				if (v== null)
					continue;
				int[] a= new int[3];
				a[0]= Integer.parseInt(ss[1]);
				a[1]= Integer.parseInt(ss[2]);
				int p= Collections.binarySearch(v, a, c);
				if (p< 0)
					p= -(p+1);
				String id= null;
				for (int i = p; i < v.size(); ++i) {
					int[] b= v.get(i);	// probe
					id= null;
					if (a[0]>b[1]|| a[1]<b[0])
						break;
					if (a[0]<= b[1]&& b[0]<=a[1]) {
						String[] start= ss[11].split(",");
						String[] len= ss[10].split(",");
						boolean stop= false;
						for (int j = 0; j < len.length&& !stop; j++) {
							int lenx= Integer.parseInt(len[j]);
							int startx= Integer.parseInt(start[j]);
							for (int k = 0; k < lenx; k++) {
								int x= a[0]+ startx+ k;
								if (x>= b[0]&& x<= b[1]) {
									if (id== null)
										id= ss[0]+":"+b[0]+"-"+b[1]; 
									String t= "";
									if (map2.containsKey(id))
										t= map2.get(id);
									map2.put(id, t+ ss[3]+ ",");
									stop= true;	// hit only once
									break;
								}
								
							}
						}
					}
				}
				
			}
			buffy.close();
			
			Iterator<String> iter= map2.keySet().iterator();
			while(iter.hasNext()) {
				String id= iter.next();
				System.out.println(id+"\t"+map2.get(id));
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
