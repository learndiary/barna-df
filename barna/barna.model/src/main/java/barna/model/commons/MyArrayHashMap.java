/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
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
