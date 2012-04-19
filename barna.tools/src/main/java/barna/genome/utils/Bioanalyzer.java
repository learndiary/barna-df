package barna.genome.utils;/*
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

import com.mindprod.ledatastream.LEDataInputStream;

import java.io.ByteArrayInputStream;
import java.io.EOFException;
import java.io.IOException;

public class Bioanalyzer {

	/**
	 * Expert Software XML Schema Date printed: 3/18/11 Agilent Technologies
	 * Page 15 of 54 Author: Volker von Einem
	 */
	static void originalSnippet() {
		
		String localName= null, content= null;
		
		if (localName.equals("RawSignal")) {
			
			ByteArrayInputStream bAIS = new ByteArrayInputStream(
					content.getBytes());
			
			Base64.InputStream b64IStream = new Base64.InputStream(bAIS,
					Base64.DECODE);
			
			LEDataInputStream lEDIStream = new LEDataInputStream(b64IStream);
			
			double yValue = 0;
			boolean endOfStream = false;
			while (!endOfStream) {
				try {
					// read values bit by bit
					yValue = (double) lEDIStream.readFloat();
					// use values
				} catch (EOFException e) {
					endOfStream = true;
				} catch (IOException e) {
					; // was uncaught in snippet
				}
			}
		}
	}
}
