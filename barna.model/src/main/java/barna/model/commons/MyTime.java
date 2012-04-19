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

package barna.model.commons;

import java.util.Calendar;

public class MyTime {

	long time= -1l;

	public MyTime(long newTime) {
		this.time= newTime;
	}
	
	public int getHours() {
		return (int) (time/ 3600000);
	}
	
	private long getRestHours() {
		return time% 3600000;
	}
	
	public int getMinutes() {
		return (int) (getRestHours()/ 60000);
	}

	private long getRestMinutes() {
		return getRestHours()% 60000;
	}
	
	public int getSeconds() {
		return (int) (getRestMinutes()/ 1000);
	}

	public long getMillis() {
		return getRestMinutes()% 1000;
	}
	
	public String toString() {
		return getHours()+ ":"+ getMinutes()+":"+ getSeconds()+"."+getMillis();
	}
	
	public static String toTime(long temps) {
		MyTime t= new MyTime(temps);
		return t.toString();
	}
	
	public static String getHexDate() {
		Calendar rightNow= Calendar.getInstance();
		StringBuffer sb= new StringBuffer(6);
		int x= rightNow.get(Calendar.YEAR);
		String s= java.lang.Integer.toString(x);
		sb.append(s.substring(2));		
		x= rightNow.get(Calendar.MONTH);
		++x;
		s= java.lang.Integer.toString(x);
		if (s.length()== 1)
			sb.append("0");
		sb.append(s);
		x= rightNow.get(Calendar.DAY_OF_MONTH);
		s= java.lang.Integer.toString(x);
		if (s.length()== 1)
			sb.append("0");
		sb.append(s);
		return sb.toString();
	}
}
