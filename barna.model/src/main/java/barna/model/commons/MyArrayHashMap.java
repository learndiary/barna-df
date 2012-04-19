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

import java.util.Arrays;

public class MyArrayHashMap<K, V> extends MyHashMap<K, V> {
	
	@Override
	public V get(Object key) {
		
			if (!(key instanceof byte[])) {
				System.err.println("only byte[] allowed");
				return null;
			}
			byte[] a= (byte[]) key;
			
	       if (key == null)
	            return getForNullKey();
	        int hash = hash(hush(a));
	        for (Entry<K,V> e = table[indexFor(hash, table.length)];
	             e != null;
	             e = e.next) {
	            Object k;
	            if (e.hash == hash && ((k = e.key) == key || iquals(a, ((byte[]) k))))	// key.equals(k))
	                return e.value;
	        }
	        return null;
	}
	
	protected int hush(byte[] a) {
	    int h= 0;
        for (int i = 0; i < a.length; i++) {
            h = 31*h + a[i];
        }
        
        return h;
	}

	/**
	 * Associates the specified value with the specified key in this map.
	 * If the map previously contained a mapping for the key, the old
	 * value is replaced.
	 *
	 * @param key key with which the specified value is to be associated
	 * @param value value to be associated with the specified key
	 * @return the previous value associated with <tt>key</tt>, or
	 *         <tt>null</tt> if there was no mapping for <tt>key</tt>.
	 *         (A <tt>null</tt> return can also indicate that the map
	 *         previously associated <tt>null</tt> with <tt>key</tt>.)
	 */
	public V put(K key, V value) {
		if (!(key instanceof byte[])) {
			System.err.println("only byte[] allowed");
			return null;
		}
		byte[] a= (byte[]) key;

	    if (key == null)
	        return putForNullKey(value);

	    int hash = hash(hush(a));
	    int i = indexFor(hash, table.length);
	    for (Entry<K,V> e = table[i]; e != null; e = e.next) {
	        Object k;
	        if (e.hash == hash && 
	        		((k = e.key) == key || 
	        				iquals(a, (byte[]) k))) {
	            V oldValue = e.value;
	            e.value = value;
	            e.recordAccess(this);
	            return oldValue;
	        }
	    }
	
	    modCount++;
	    addEntry(hash, key, value, i);
	    return null;
	}
	
	boolean iquals(byte[] a, byte[] b) {
		return Arrays.equals(a, b);
	}
	

	/**
	 * Constructs an empty <tt>HashMap</tt> with the specified initial
	 * capacity and the default load factor (0.75).
	 *
	 * @param  initialCapacity the initial capacity.
	 * @throws IllegalArgumentException if the initial capacity is negative.
	 */
	public MyArrayHashMap(int initialCapacity) {
	    super(initialCapacity);
	}
	
	
}
