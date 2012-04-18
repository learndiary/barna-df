import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;


public class RgaspCracker {
	public static void main(String[] args) {
		String fName= "N:\\rgasp1.2\\postrgasp\\stat3_chr17_40465342-40540513C.txt";
		try {
			System.err.println("LOCUS: "+ fName);
			System.err.println("======\n");
			BufferedReader buffy= new BufferedReader(new FileReader(fName));
			String s;
			HashMap<Float, String> mapExpr= new HashMap<Float, String>();
			while ((s= buffy.readLine())!= null) {
				String[] ss= s.split("\\s");
				System.err.println(ss[0]+"\t"+ss[1]);
				Float v= Float.parseFloat(ss[1]);
				String id= mapExpr.get(v);
				if (id== null)
					id= "";
				else
					id+= "/";
				id+= ss[0]; 
				mapExpr.put(v, id);
			}
			buffy.close();
			System.err.println();
			
			crack(mapExpr, 23.9f);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	static strictfp void crack(HashMap<Float, String> mapExpr, float val) {
		
		Float zero= new Float(0.0);
	
		// make non-zero spectrum
		int dim= mapExpr.size();
		if (mapExpr.containsKey(zero))
			--dim;
		float[] spec= new float[dim];
		Iterator<Float> iter= mapExpr.keySet().iterator();
		int p= 0;
		while (iter.hasNext()) {
			float x= iter.next();
			if (x== 0)
				continue;
			spec[p++]= x;
		}
		
		// crack
		boolean[] acc= new boolean[spec.length];
		for (int i = 0; i < acc.length; i++) 
			acc[i]= false;
		HashMap<Float, Vector<boolean[]>> mapRes= new HashMap<Float, Vector<boolean[]>>();
		for (int i = 0; i < spec.length; i++) {
			crackRek(spec, i, val, 0, acc.clone(), mapRes);
		}
		
		// get best results
		System.err.println("BEST SOLUTION: "+val);
		System.err.println("==============\n");
		float[] dd= new float[mapRes.size()];
		iter= mapRes.keySet().iterator();
		p= 0;
		while (iter.hasNext())  
			dd[p++]= iter.next();
		Arrays.sort(dd);
		for (int i = 0; i < Math.min(4,dd.length); i++) {
			if (false&& i> 0&& (Math.abs(dd[i]- val)/ val)> 1)
				break;
			Vector<boolean[]> v= mapRes.get(dd[i]);
			for (int x = 0; x < v.size(); ++x) {
				System.err.print(dd[i]+"\t");
				StringBuilder sb= new StringBuilder();
				boolean[] bb= v.elementAt(x);
				float sum= 0;
				for (int j = 0; j < bb.length; j++) {
					if (!bb[j]) 
						continue;
					sb.append(mapExpr.get(spec[j])+ ",");
					sum+= spec[j];
					System.err.print(spec[j]+ "+");
				}
				System.err.println(" = "+ sum+ "("+ Math.round(sum*100)/ 100f+ ")");
				System.err.println(sb);
			}
		}
		if (mapExpr.containsKey(zero)) {
			System.err.println("\nNULL:");
			System.err.println("=====");
			System.err.println(mapExpr.get(zero).replaceAll("/", ","));
		}
			
	}
	
	static void crackRek(float[] spec, int p, float target, float sum, boolean[] a, HashMap<Float,Vector<boolean[]>> mapRes) {
		
		if (Math.round((sum+ spec[p])* 100)/100f>= target|| p+1== spec.length) {
			float dlo= Math.abs((Math.round(sum* 100)/100f)- target),
					dhi= Math.abs((Math.round((sum+ spec[p])* 100)/100f)- target);
			
			float min= Math.min(dlo,dhi);
			if (min== dhi)
				a[p]= true;
			Vector<boolean[]> v= mapRes.get(min);
			if (v== null)
				v= new Vector<boolean[]>();
			v.add(a);
			mapRes.put(min,v);
			return;
		}
		
		a[p]= true;
		sum+= spec[p];
		for (int i = p+1; i < spec.length; i++) {
			crackRek(spec, i, target, sum, a.clone(), mapRes);
		}
		
	}
}
