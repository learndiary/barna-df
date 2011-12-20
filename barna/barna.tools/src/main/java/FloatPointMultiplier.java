import barna.genome.model.commons.MyFile;

import java.io.*;


public class FloatPointMultiplier {

	public static void main(String[] args) {
		replace(new File("P:\\rgasp1.3\\GM12878\\quant01_fixed-25lines_others_dmg\\GM12878_flux_paired.gtf"), -1, 1.48f);
	}

	private static void replace(File file, int nr, float f) {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(file));
			String oFname= MyFile.stripExtension(file.getAbsolutePath())+ "_x_"+ Float.toString(f);
			BufferedWriter writer= new BufferedWriter(new FileWriter(oFname));
			String line;
			int cnt= 0;
			while ((line= buffy.readLine())!= null) {
				++cnt;
				String[] ss= line.split("\\s");
				String token= ss[ss.length- 1];
				float b= Float.parseFloat(token.substring(0, token.length()- 1));
				b*= f;
				line= line.substring(0, line.length()- ss[ss.length- 1].length());
				writer.write(line+ String.format("%1$f", b)+";\n");
			}
			buffy.close();
			writer.flush();
			writer.close();
			System.err.println("wrote "+ cnt+ " lines.");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
