package barna.astafunk.DP;

/**
 * @version 2.0
 * @author Vitor Coelho
 * @since Last release in 2014-07-17
 */

/**
 * This class describes a cell of a dynamic programming matrix.
 */
public class DyPCell {

    /**
     * The previous cell of the track.
     */
    private DyPCell prevCell;

    /**
     * Bit score of the DP matrix cell.
     */
    private double score = 0;

    /**
     * Score matrix row.
     */
    private int row;

    /**
     * Score matrix column.
     */
    private int col;

    /**
     * Constructor of Dynamic Programming matrix cell.
     * @param row Number of row
     * @param col Number of column
     */
    public DyPCell(int row, int col) {

        this.col = col;
        this.row = row;
    }

    /**
     *
     * @return Previous cell.
     */
    public DyPCell getPrevCell() {
        return prevCell;
    }

    /**
     *
     * @param prevCell Previous cell to set.
     */
    public void setPrevCell(DyPCell prevCell) {
        this.prevCell = prevCell;
    }

    /**
     *
     * @return Bit score of the cell.
     */
    public double getScore() {
        return score;
    }

    /**
     * M
     * @param score Bit score (log 2) of the cell to set.
     */
    public void setScore(double score) {
        this.score = score;
    }

    /**
     *
     * @return Row number of the cell.
     */
    public int getRow() {
        return row;
    }

    /**
     *
     * @return Column number of the cell.
     */
    public int getCol() {
        return col;
    }

}

