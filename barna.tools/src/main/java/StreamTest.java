
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

public class StreamTest {

	public static void main(String[] args) {
		for (int i = 0; i < 10; i++) {
			System.err.print("working [");
			for (int j = 0; j < i; j++) 
				System.err.print('-');
			System.err.print('\\');
			for (int j = i+1; j < 10; j++) 
				System.err.print(' ');
			System.err.print("]\r");
			System.err.flush();
			
			try {
				Thread.currentThread().sleep(100);
			} catch (InterruptedException e) {
				; // :)
			}

			System.err.print("working [");
			for (int j = 0; j < i; j++)
				System.err.print('-');
			System.err.print('|');
			for (int j = i+1; j < 10; j++) 
				System.err.print(' ');
			System.err.print("]\r");
			System.err.flush();
			
			try {
				Thread.currentThread().sleep(100);
			} catch (InterruptedException e) {
				; // :)
			}

			System.err.print("working [");
			for (int j = 0; j < i; j++)
				System.err.print('-');
			System.err.print('/');
			for (int j = i+1; j < 10; j++) 
				System.err.print(' ');
			System.err.print("]\r");
			System.err.flush();
			
			try {
				Thread.currentThread().sleep(100);
			} catch (InterruptedException e) {
				; // :)
			}

			System.err.print("working [");
			for (int j = 0; j <= i; j++) 
				System.err.print('-');
			for (int j = i+1; j < 10; j++) 
				System.err.print(' ');
			System.err.print("]\r");
			System.err.flush();
			
			try {
				Thread.currentThread().sleep(100);
			} catch (InterruptedException e) {
				; // :)
			}
			
		}
	}
}
