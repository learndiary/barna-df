/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

import barna.model.commons.MyFile;

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
