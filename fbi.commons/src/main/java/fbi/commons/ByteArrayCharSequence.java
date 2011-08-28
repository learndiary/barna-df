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

package fbi.commons;

import java.util.Arrays;

/**
 * Mutable string implementation with field support
 *
 * @author Micha Sammeth (gmicha@googlemail.com)
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ByteArrayCharSequence implements CharSequence, Comparable<CharSequence> {
    /**
     * Tab character byte value
     */
    public static final byte BYTE_TAB = 9;
    /**
     * DNA/RNA characters
     */
    private static final byte[] CHARS_NORMAL = new byte[]{
            45, 65, 67, 71, 75, 77, 78, 82, 84, 85, 88, 89,        // -ACGKMNRTUXY
            97, 99, 103, 107, 109, 110, 114, 116, 117, 120, 121};// acgkmnrtuxy
    /**
     * DNA/RNA complements
     */
    private static final byte[] CHARS_REVERSED = new byte[]{
            45, 84, 71, 67, 77, 75, 78, 89, 65, 65, 88, 82,        // -TGCMKNYAAXR
            116, 103, 99, 109, 107, 110, 121, 97, 97, 120, 114};    // tgcmknyaaxr


    private final static byte[] DigitTens = {
            48, 48, 48, 48, 48, 48, 48, 48, 48, 48,
            49, 49, 49, 49, 49, 49, 49, 49, 49, 49,
            50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
            51, 51, 51, 51, 51, 51, 51, 51, 51, 51,
            52, 52, 52, 52, 52, 52, 52, 52, 52, 52,
            53, 53, 53, 53, 53, 53, 53, 53, 53, 53,
            54, 54, 54, 54, 54, 54, 54, 54, 54, 54,
            55, 55, 55, 55, 55, 55, 55, 55, 55, 55,
            56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
            57, 57, 57, 57, 57, 57, 57, 57, 57, 57,
    };

    private final static byte[] DigitOnes = {
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
    };


    final static int[] sizeTable = {9, 99, 999, 9999, 99999, 999999, 9999999,
            99999999, 999999999, Integer.MAX_VALUE};

    /**
     * All possible chars for representing a number as a String
     */
    private final static byte[] digits = {
            48, 49, 50, 51, 52, 53,
            54, 55, 56, 57, 97, 98,
            99, 100, 101, 102, 103, 104,
            105, 106, 107, 108, 109, 110,
            111, 112, 113, 114, 115, 116,
            117, 118, 119, 120, 121, 122
    };
    /**
     * The characters array
     */
    public byte[] chars;
    /**
     * Start index
     */
    public int start;
    /**
     * End index (excluded, first index NOT in the string)
     */
    public int end;    // end excluded, heilichs blechle !!
    /**
     * Helper for internal find
     */
    protected int p1;
    /**
     * Helper for internal find
     */
    protected int p2;
    /**
     * Helper for internal find
     */
    protected int cnt;

    /**
     * Separator character
     */
    private char separator = '\t';
    /**
     * The load factor used to excend the buffer size
     */
    private double loadFactor = 0.3;

    /**
     * Initialize the instance using the given byte buffer. The buffers reference is
     * used directly, no copy is created and start and end ar set to 0 and length of the
     * buffer respectively
     *
     * @param value the byte buffer
     */
    public ByteArrayCharSequence(byte[] value) {
        if(value == null) throw new NullPointerException("Null buffer not allowed!");
        chars = value;
        this.start = 0;
        this.end = chars.length;
        resetFind();
    }

    /**
     * Create an instance from the given string
     *
     * @param s the string
     */
    public ByteArrayCharSequence(String s) {
        this(s.getBytes());
    }

    /**
     * Create a new empty instance with the initial capacity
     *
     * @param capacity the initial capacity
     */
    public ByteArrayCharSequence(int capacity) {
        this(new byte[capacity]);
        this.end = 0;
        resetFind();
    }

    /**
     * Clone/Copy constructor to create a clone of the given sequence
     *
     * @param source the source
     */
    public ByteArrayCharSequence(ByteArrayCharSequence source) {
        this(source.cloneBuffer());
        /*
        Clone buffer reduces the chars array already to the string bounded by start and end
        so we have to use 0 and length as new bounds
         */
        this.start = 0;
        this.end = chars.length;
        this.separator = source.separator;

    }

    /**
     * INTERNAL: creates a new instance with given bounds
     *
     * @param value the values
     * @param start the start inclusive
     * @param end the end exclusive
     */
    protected ByteArrayCharSequence(byte[] value, int start, int end) {
        this(value);
        this.start = start;
        this.end = end;
    }


    /**
     * Set the value to the given string
     *
     * @param value the string
     */
    public void set(String value) {
        if(value == null) throw new NullPointerException("Null value not permitted!");
        chars = value.getBytes();
        this.start = 0;
        this.end = chars.length;
        resetFind();
    }

    /**
     * Get char value at given index
     *
     * @param index the index
     * @return character the character at the specivied index
     */
    @Override
    public char charAt(int index) {
        return (char) chars[start + index];
    }

    /**
     * Get byte value of the character at the specified index
     *
     * @param index the index
     * @return byte byte value at the specified index
     */
    public byte byteAt(int index) {
        return chars[start + index];
    }

    /**
     * Returns the length of the current sequence
     *
     * @return length the length of the current sequence
     */
    @Override
    public int length() {
        return (end - start);
    }


    /**
     * Remove all white-spaces from start and end
     */
    public void trim() {
        while ((start < end) && (chars[start] <= ' ')) {
            start++;
        }
        while ((start < end) && (chars[end - 1] <= ' ')) {
            end--;
        }
    }

    /**
     * Create a sub-sequence of the current sequence. Note that
     * the sub-sequence operates on the same byte buffer! Changes to the
     * buffer affect this instance and the sub-sequence
     *
     * @param start start index inclusive
     * @param end end index exclusive
     * @return subsequence the subsequence
     */
    @Override
    public ByteArrayCharSequence subSequence(int start, int end) {
        return new ByteArrayCharSequence(chars, this.start + start, this.start
                + end);
    }


    /**
     * Search for a given field. Fields are separated by the {@link #separator}.
     * Returns true if the field was found. The field variables {@link #p1} and {@link #p2}
     * are set to the left/right bound of the field.
     *
     *
     * @param fieldNr 0-based field index
     * @return found return true if the field was found
     */
    protected boolean find(int fieldNr) {
        if(fieldNr < 0) return false;
        if (fieldNr == cnt) {
            return (p1 != p2 || p1 != end);
        }

        if (fieldNr > cnt) {
            if (p1 != end) {
                for (int i = p2 + 1; cnt < fieldNr && i < end; ++i) {
                    if (charAt(i) == separator) {
                        ++cnt;
                        p1 = p2 + 1;
                        p2 = i;
                    }
                }
                //++p1;
                if (cnt < fieldNr) {
                    p1 = (p2 == end) ? p2 : p2 + 1;
                    p2 = end;
                    ++cnt;
                }
            }
        } else {
            if (p1 != start) {
                for (int i = p1 - 2; cnt > fieldNr && i >= start; --i) {
                    if (charAt(i) == separator) {
                        --cnt;
                        p2 = p1 - 1;
                        p1 = i + 1;
                    }
                }
                if (cnt > fieldNr) {
                    p2 = (p1 == start) ? start : p1 - 1;
                    p1 = start;
                    --cnt;
                }
            }
        }
        return cnt == fieldNr && (p1 != p2 || p1 != end);
    }

    /**
     * Get the token at the given field. This returns null if the field is
     * not found, otherwise {@link #subSequence(int, int)}  is used to
     * create the token representation. Note that {@link #subSequence(int, int)} creates
     * a view on the sequence, therefore the same byte buffer is used and changes
     * affect both this instance and the returned token
     *
     * @param fieldNr the field
     * @return token the token or null
     */
    public ByteArrayCharSequence getToken(int fieldNr) {
        if (fieldNr < 0) {
            return null;
        }
        if(!find(fieldNr)) return null;
        return subSequence(p1, p2);
    }

    /**
     * Get the field at the given index and parse it to an Integer.
     * An exception is triggered if either the field does not exist or
     * the value could not be parsed
     *
     * @param fieldNr the field index
     * @return value integer value of the field
     */
    public int getTokenInt(int fieldNr) {
        if (fieldNr < 0) {
            throw new IllegalArgumentException("Field index < 0");
        }

        if(find(fieldNr)){
            return parseInt(p1, p2);
        }
        throw new IllegalArgumentException("Field " + fieldNr + " not found");
    }

    /**
     * Replace the value of the given field with the given integer
     *
     * @param fieldNr the field index
     * @param value the value
     * @return success true if replaced
     */
    public boolean replace(int fieldNr, int value) {
        return fieldNr >= 0 && find(fieldNr) && replaceCurrField(value);
    }

    /**
     * Replace the given field with the give character sequence
     *
     * @param fieldNr the field index
     * @param value the value
     * @return success true if success
     */
    public boolean replace(int fieldNr, CharSequence value) {
        return fieldNr >= 0 && find(fieldNr) && replaceCurrField(value);

    }

    /**
     * Replace the current field with the given byte
     *
     * @param value the value
     * @return true if success
     */
    protected boolean replaceCurrField(byte value) {
        if (p1 < start || p1 >= end) {
            return false;
        }

        int diff = 1 - (p2 - p1);
        if (end + diff > chars.length) {
            extend(diff);
        }
        if (diff != 0) {
            System.arraycopy(chars, p2, chars, p2 + diff, end - p2);
            end += diff;
            p2 += diff;
        }
        assert (p2 == p1 + 1);
        chars[p1] = value;
        return true;
    }

    /**
     * Replace the current field (bounds stored in {@link #p1} and {@link #p2}) with the given
     * integer value
     *
     * @param value the value
     * @return success true if replaced
     */
    protected boolean replaceCurrField(int value) {
        if (p1 < start || p1 > end){
            return false;
        }

        int digits = value < 0 ? countCharacters(-value) + 1 : countCharacters(value);
        int diff = digits - (p2 - p1);
        if (end + diff > chars.length) {
            extend(diff);
        }
        if (diff != 0) {
            System.arraycopy(chars, p2, chars, p2 + diff, end - p2);
            end += diff;
            p2 += diff;
        }

        int p = insertAsCharacters(value, p2, chars);    // p2 is new end of field
        assert (p == p1);
        return true;
    }

    /**
     * Replace the current field (bounds stored in {@link #p1} and {@link #p2}) with the given
     * character sequence
     *
     * @param value the value
     * @return success true if replaced
     */
    protected boolean replaceCurrField(CharSequence value) {
        if (p1 < start || p1 > end)    // for search ends at start/end, p1== p2
        {
            return false;
        }

        int digits = value.length();
        int diff = digits - (p2 - p1);
        if (end + diff > chars.length) {
            extend(diff);
        }
        if (diff != 0) {
            System.arraycopy(chars, p2, chars, p2 + diff, end - p2);
            end += diff;
            p2 += diff;
        }

        int p;
        for (p = p1; p < p1 + digits; ++p) {
            chars[p] = (byte) value.charAt(p - p1);
        }

        assert (p == p2);
        return true;
    }

    @Override
    public int hashCode() {
        int h = 0, len = length();

        for (int i = 0; i < len; i++) {
            h = 31 * h + charAt(i);
        }

        return h;
    }

    /**
     * Append the given character sequence to the end
     *
     * @param cs the characters to append
     */
    public void append(CharSequence cs) {
        int len = cs.length();
        ensureLength(end, len);
        for (int i = 0; i < len; ++i) {
            chars[end + i] = (byte) cs.charAt(i);
        }
        end += len;
    }

    /**
     * Append the content of the given byte array limited by from and to indices
     * to this sequence
     *
     * @param b the bytes to be appended
     * @param from start index inclusive
     * @param to end index exclusive
     */
    public void append(byte[] b, int from, int to) {
        int len = to - from;
        ensureLength(end, len);
        System.arraycopy(b, from, chars, end, len);
        end += len;
    }

    /**
     * Append the given byte to this sequence
     *
     * @param x the byte to be appended
     */
    public void append(byte x) {
        ensureLength(end, 1);
        chars[end++] = x;
    }

    /**
     * Append the given char to this sequence
     * @param x the char
     */
    public void append(char x) {
        ensureLength(end, 1);
        chars[end++] = (byte)x;
    }

    /**
     * Append the given integer to this sequence
     *
     * @param x the integer
     */
    public void append(int x) {
        int len = countCharacters(x);
        ensureLength(end, len);
        end += len;    // end index
        insertAsCharacters(x, end, chars);
    }

    /**
     * Append the given byte to the current field
     *
     * @param x the byte to be appended
     */
    protected void appendCurrField(byte x) {
        if(!find(cnt)) throw new IllegalArgumentException("Current field " +cnt + " is not valid!");
        ensureLength(end, 1);
        System.arraycopy(chars, p2, chars, p2 + 1, end - p2);
        end++;
        chars[p2++] = x;
    }

    /**
     * Append the given integer to the current field
     *
     * @param x the integer
     */
    protected void appendCurrField(int x) {
        if(!find(cnt)){
            //throw new IllegalArgumentException("Current field " +cnt + " is not valid!");

        }
        int digits = countCharacters(x);
        ensureLength(end, digits);
        System.arraycopy(chars, p2, chars, p2 + digits, end - p2);
        end += digits;
        p2 += digits; // end index
        insertAsCharacters(x, p2, chars);
    }

    /**
     * Append the given sequence to the current field
     * @param value the sequence
     */
    protected void appendCurrField(CharSequence value) {
        if(value == null) throw new NullPointerException();
        if(!find(cnt)) throw new IllegalArgumentException("Current field " +cnt + " is not valid!");
        int digits = value.length();
        ensureLength(end, digits);
        if (digits != 0) {
            System.arraycopy(chars, p2, chars, p2 + digits, end - p2);
            end += digits;
            p2 += digits;
        }

        int p;
        for (p = p1; p < p1 + digits; ++p) {
            chars[p] = (byte) value.charAt(p - p1);
        }

        assert (p == p2);
    }

    /**
     * Extend the buffer size based on the load factor
     */
    // todo: this should be protected and not used from outside
    public void extend() {
        extend(Math.max(1, (int) (chars.length * loadFactor)));
    }

    /**
     * Extend the buffer size
     *
     * @param x extend
     */
    public void extend(int x) {
        byte[] b = new byte[chars.length + x];
        System.arraycopy(chars, 0, b, 0, chars.length);
        chars = b;
    }

    /**
     * Creates a copy of this sequence
     *
     * @return copy the copy
     */
    public ByteArrayCharSequence cloneCurrentSeq() {
        return new ByteArrayCharSequence(this);
    }

    /**
     * Returns a copy of the current buffer limited to the string representation
     * defined by {@link #start} and {@link #end}.
     *
     * @return clone clone of the buffer limited to {@link #start} and {@link #end}.
     */
    protected byte[] cloneBuffer() {
        byte[] cloned = new byte[length()];
        System.arraycopy(chars, start, cloned, 0, cloned.length);
        return cloned;
    }

    @Override
    public String toString() {
        return new String(chars, start, length());
    }

    /**
     * Returns a char array representing this sequence
     *
     * @return chars the char array representing this sequence
     */
    public char[] toCharArray() {
        char[] cc = new char[length()];
        for (int i = 0; i < cc.length; i++) {
            cc[i] = charAt(i);
        }
        return cc;
    }

    /**
     * Copies the current sequence content to the given char array. If the
     * given array is not large enough, a new array is created
     *
     * @param c the target (only used if large enough)
     * @return chars the char array representing the current sequence
     */
    public char[] toCharArray(char[] c) {
        if (length() > c.length) {
            return toCharArray();
        }
        for (int i = 0; i < length(); i++) {
            c[i] = charAt(i);
        }
        return c;
    }

    /**
     * Returns true if this sequence starts with the given sequence
     *
     * @param cs the sequence
     * @return startsWith true if this sequence starts with the given sequence
     */
    public boolean startsWith(CharSequence cs) {
        if(cs == null) throw new NullPointerException();
        int length = cs.length();
        if (length() < length) {
            return false;
        }
        for (int i = 0; i < length; i++) {
            if (cs.charAt(i) != charAt(i)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns true if this sequence ends with the given sequence
     *
     * @param cs the sequence
     * @return startsWith true if this sequence starts with the given sequence
     */
    public boolean endsWith(CharSequence cs) {
        if(cs == null) throw new NullPointerException();
        int length = cs.length();
        int thisLength = length();
        if (thisLength < length) {
            return false;
        }
        for (int i = 0; i < length; i++) {
            if (cs.charAt(length - i-1) != charAt(thisLength - i-1)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Parse the integer value for this sequence
     *
     * @return value the integer value
     */
    public int parseInt() {
        int val = 0;
        for (int i = end - 1, pow = 1; i >= start; --i, pow *= 10) {
            if(Arrays.binarySearch(digits,chars[i] ) < 0){
                // not a digit
                throw new NumberFormatException("Unable to parse " + toString() + " to an integer!");
            }
            val += (chars[i] - 48) * pow;
        }
        return val;
    }

    /**
     * Parse the integer value for the sequence between
     * from (inclusive) and to (exclusive)
     *
     * @param from the start index (inclusive)
     * @param to the end index (exclusive)
     * @return value the integer value
     */
    public int parseInt(final int from, final int to) {
        if(from == to) throw new NumberFormatException("Unable to parse an empty string to a number");
        int f = 1;
        int start = from;

        // check sign
        if (chars[from] == 45) {
            f = -1;
            ++start;
        }


        int val = 0;
        for (int i = to - 1, pow = 1; i >= start; --i, pow *= 10) {
            // check if its a digit
            if(Arrays.binarySearch(digits,chars[i] ) < 0){
                // not a digit
                throw new NumberFormatException("Unable to parse " + subSequence(from, to) + " to an integer!");
            }
            val += (chars[i] - 48) * pow;
        }
        return (val * f);
    }

    /**
     * Compares two strings lexicographically.
     * The comparison is based on the Unicode value of each character in
     * the strings. The character sequence represented by this
     * <code>String</code> object is compared lexicographically to the
     * character sequence represented by the argument string. The result is
     * a negative integer if this <code>String</code> object
     * lexicographically precedes the argument string. The result is a
     * positive integer if this <code>String</code> object lexicographically
     * follows the argument string. The result is zero if the strings
     * are equal; <code>compareTo</code> returns <code>0</code> exactly when
     * the {@link #equals(Object)} method would return <code>true</code>.
     * <p/>
     * This is the definition of lexicographic ordering. If two strings are
     * different, then either they have different characters at some index
     * that is a valid index for both strings, or their lengths are different,
     * or both. If they have different characters at one or more index
     * positions, let <i>k</i> be the smallest such index; then the string
     * whose character at position <i>k</i> has the smaller value, as
     * determined by using the &lt; operator, lexicographically precedes the
     * other string. In this case, <code>compareTo</code> returns the
     * difference of the two character values at position <code>k</code> in
     * the two string -- that is, the value:
     * <blockquote><pre>
     * this.charAt(k)-anotherString.charAt(k)
     * </pre></blockquote>
     * If there is no index position at which they differ, then the shorter
     * string lexicographically precedes the longer string. In this case,
     * <code>compareTo</code> returns the difference of the lengths of the
     * strings -- that is, the value:
     * <blockquote><pre>
     * this.length()-anotherString.length()
     * </pre></blockquote>
     *
     * @param anotherSeq the <code>String</code> to be compared.
     * @return the value <code>0</code> if the argument string is equal to
     *         this string; a value less than <code>0</code> if this string
     *         is lexicographically less than the string argument; and a
     *         value greater than <code>0</code> if this string is
     *         lexicographically greater than the string argument.
     */
    @Override
    public int compareTo(CharSequence anotherSeq) {
        int len1 = length();
        int len2 = anotherSeq.length();
        int n = Math.min(len1, len2);

        for (int m = 0; m < n; m++) {
            if (charAt(m) != anotherSeq.charAt(m)) {
                return charAt(m) - anotherSeq.charAt(m);
            }
        }
        return len1 - len2;

    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if(obj instanceof CharSequence){
            CharSequence cs = (CharSequence) obj;
            if (cs.length() != length()) {
                return false;
            }
            for (int i = 0; i < cs.length(); i++) {
                if (cs.charAt(i) != charAt(i)) {
                    return false;
                }
            }
            return true;
        }
        return false;
    }

    /**
     * Count the number of available tokens
     *
     * @return tokens the number of tokens
     */
    public int countTokens() {
        int cnt = 0;
        for (int i = start; i < length(); i++) {
            if (charAt(i) == separator) {
                ++cnt;
            }
        }
        if (length() > 0 && charAt(length() - 1) != separator) {
            ++cnt;
        }
        return cnt;
    }

    /**
     * Reset find caches
     */
    public void resetFind() {
        p1 = start - 1;
        p2 = start - 1;
        cnt = -1;
    }

    /**
     * Clear this sequence.
     */
    public void clear() {
        start = end = 0;
        resetFind();
    }

    /**
     * Ensure that the buffer can contain data from the given index with
     * given length.
     *
     * @param from the index
     * @param len the length
     */
    public void ensureLength(int from, int len) {
        int diff = (from + len) - chars.length;
        if (diff > 0) {
            int e = chars.length * 2;
            if (e > diff) {
                diff = e;
            }
            extend(diff);
        }
    }

    /**
     * Convert the sequence to upper case
     *
     */
    public void toUpperCase() {
        toUpperCase(chars, start, end);
    }
    /**
     * Convert the character between {@code from}(inclusive) and {@code to}(exclusive) to upper case
     *
     * @param from start index (inclusive)
     * @param to end index (exclusive)
     */
    public void toUpperCase(int from, int to) {
        toUpperCase(chars, start + from, start + to);
    }
    /**
     * Convert the sequence to lower case
     *
     */
    public void toLowerCase() {
        toLowerCase(chars, start, end);
    }
    /**
     * Convert the character between {@code from}(inclusive) and {@code to}(exclusive) to lower case
     *
     * @param from start index (inclusive)
     * @param to end index (exclusive)
     */
    public void toLowerCase(int from, int to) {
        toLowerCase(chars, start + from, start + to);
    }

    /**
     * Convert the given sequence to its reverse complement
     *
     * @param cs the sequence
     * @param from start index (inclusive)
     * @param to end index (exclusive)
     */
    public static void reverseComplement(ByteArrayCharSequence cs, int from, int to) {
        complement(cs.chars, cs.start + from, cs.start + to);
        reverse(cs.chars, cs.start + from, cs.start + to);
    }
    /**
     * Convert the given sequence to its reverse complement
     *
     * @param cs the sequence
     * @param from start index (inclusive)
     * @param to end index (exclusive)
     */
    public static void reverse(ByteArrayCharSequence cs, int from, int to) {
        reverse(cs.chars, cs.start + from, cs.start + to);
    }

    /**
     * Reverse the content of the byte array between {@code from} and {@code to-1}
     *
     * @param a the array
     * @param from the start index (inclusive)
     * @param to the end index (exclusive)
     */
    public static void reverse(byte[] a, int from, int to) {
        // adapted from commons.ArrayUtils
        if (a == null) {
            return;
        }
        int i = from;
        int j = to - 1;
        byte tmp;
        while (j > i) {
            tmp = a[j];
            a[j] = a[i];
            a[i] = tmp;
            --j;
            ++i;
        }

    }

    /**
     * Convert the given array to its complement
     *
     * @param a the array
     * @param from start index (inclusive)
     * @param to end index (exclusive)
     */
    public static void complement(byte[] a, int from, int to) {
        int p;
        for (int i = from; i < to; ++i) {
            p = Arrays.binarySearch(CHARS_NORMAL, a[i]);
            if (p < 0) {
                throw new RuntimeException("Complement: unknown symbol " + ((char) a[i]) + " (" + a[i] + ")");
            }
            a[i] = CHARS_REVERSED[p];
        }
    }



    /**
     * Convert the character between {@code from}(inclusive) and {@code to}(exclusive) to upper case
     *
     * @param b the array
     * @param from start index (inclusive)
     * @param to end index (exclusive)
     */
    public static void toUpperCase(byte[] b, int from, int to) {
        for (int i = from; i < to; i++) {
            if (b[i] >= 97 && b[i] <= 122)    // 'a'..'z'
            {
                b[i] -= 32;
            }
        }
    }

    /**
     * Convert the character between {@code from}(inclusive) and {@code to}(exclusive) to lower case
     *
     * @param b the array
     * @param from start index (inclusive)
     * @param to end index (exclusive)
     */
    public static void toLowerCase(byte[] b, int from, int to) {
        for (int i = from; i < to; i++) {
            if (b[i] >= 65 && b[i] <= 90)    // 'A'..'Z'
            {
                b[i] += 32;
            }
        }
    }

    /**
     * Returns the number of characters needed to store the given integer in
     * its string representation.
     *
     * @param number the integer to store
     * @return characters number of characters needed to convert the given number to a string
     */
    protected static int countCharacters(int number) {
        int val = 0;
        if (number < 0) {
            val = 1;
            number = -number;
        }
        for (int i = 0; ; i++) {
            if (number <= sizeTable[i]) {
                return val + i + 1;
            }
        }
    }

    /**
     * Convert the given integer to characters and set in in the given byte array at
     * the given index
     *
     *
     * @param i the integer
     * @param endIndex the end index
     * @param buf the buffer
     * @return start start index of the integer after insert
     */
    protected static int insertAsCharacters(int i, int endIndex, byte[] buf) {
        int q, r;
        int charPos = endIndex;
        byte sign = 0;

        if (i < 0) {
            sign = 45;    // '-'
            i = -i;
        }

        // Generate two digits per iteration
        while (i >= 65536) {
            q = i / 100;
            // really: r = i - (q * 100);
            r = i - ((q << 6) + (q << 5) + (q << 2));
            i = q;
            buf[--charPos] = DigitOnes[r];
            buf[--charPos] = DigitTens[r];
        }

        // Fall thru to fast mode for smaller numbers
        // assert(i <= 65536, i);
        for (; ; ) {
            q = (i * 52429) >>> (16 + 3);
            r = i - ((q << 3) + (q << 1));  // r = i-(q*10) ...
            buf[--charPos] = digits[r];
            i = q;
            if (i == 0) {
                break;
            }
        }
        if (sign != 0) {
            buf[--charPos] = sign;
        }

        return charPos;
    }

    public static ByteArrayCharSequence cloneSequence(ByteArrayCharSequence src, ByteArrayCharSequence target) {
        if (target == null) {
            target = new ByteArrayCharSequence(src.length());
        }
        if (target.chars.length < src.length()) {
            target.chars = new byte[src.length()];
        }
        System.arraycopy(src.chars, src.start, target.chars, 0, src.length());
        target.start = src.start;
        target.end = src.end;
        return target;
    }


}
