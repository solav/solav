function waterman_eggert(gapopen_pen, gapextend_pen, seqh, seqv, scoring_matrix, threshold) 
{   
    //setting up scoring matrix and amino acid lookup
    var amino_acid_dict = {}
    var amino_acids = "A   R   N   D   C   Q  E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X".split(/\W+/)
    var i,j;
    for (i=1; i < amino_acids.length; i++) {
        amino_acid_dict[amino_acids[i-1]] = i
    }
    var score_matrix = get_scoring_matrix(scoring_matrix)

    //directions enum
    var dir = Object.freeze({NONE: 0, UP:1, LEFT:2, DIAGONAL:3}) // alignment direction
    
    //Prepend character to make string indexes line up with matrix indexes
    passing_scores = {}
    seqv = " " + seqv
    seqh = " " + seqh 
    if (gapopen_pen > 0)
        gapopen_pen = -gapopen_pen
    if (gapextend_pen > 0)
        gapextend_pen = -gapextend_pen
    //initialize matrix
    var matrix = new Array(seqv.length);
    
    for (i=0; i < matrix.length; i++) {
        matrix[i] = new Array(seqh.length);
    }

    for (i=0; i < matrix.length; i++) {
        for (j=0; j < matrix[0].length; j++) {
            matrix[i][j] = { direction: dir.NONE, score: 0 }
        }
    }

    var used_set = {}
    var top_scores = {}
    
    
    var first_pass = true
    //the arguments are where to recalculate from i.e. the start of the last solution
    function smith_waterman(srow,scol) {  
        var up, left, diagonal, up_score, left_score, diagonal_score, max_score, max_coords
        var i, previous_row_last_changed, this_row_last_changed
        
        // Limiting of recalculation will work like this:
        //  1. at the beginning set previous_row_last_changed to 1
        //  2. in current row i: for each column j 
        //      set this_row_last_changed = 0
        //      if the recalculation changes the value, set this_row_last_changed to j and continue to next column, otherwise
        //              if j-1 < previous_row_last_changed, continue to next j, (this cell won't change the next column's value but the up or diagonal direction might)
        //              else previous_row_last_changed = this_row_last_changed (no cells to the right need recalculation)
        //                  continue to the next row i (outerfor)
        //      if previous_row_last_changed == 0,  break outerloop. (We know the rest of the cells need no more recalculation.)


        previous_row_last_changed = matrix[0].length


        
        previous_row_last_changed = -1; //there were no changes in the previous row but don't trigger the break
        outerloop:
        for (i=srow; i < matrix.length; i++) {
            if (previous_row_last_changed==0) 
                break
            this_row_last_changed = 0;            
            innerloop:
            for (j=scol; j < matrix[0].length; j++) {
                        
                up = matrix[i-1][j]
                left = matrix[i][j-1]
                diagonal = matrix[i-1][j-1]
                
                if (left.direction == dir.LEFT)
                    left_score = left.score + gapextend_pen
                else 
                    left_score = left.score + gapopen_pen
 
                if (up.direction == dir.UP)
                    up_score = up.score + gapextend_pen
                else
                    up_score = up.score + gapopen_pen
 
                if ([i,j].toString() in used_set)
                    diagonal_score = -1           // remove from contention
                else
                    diagonal_score = diagonal.score + score_matrix[amino_acid_dict[seqv[i]]][amino_acid_dict[seqh[j]]]
 
                switch (Math.max(diagonal_score,0,up_score,left_score)){
                    case 0:
                        if (matrix[i][j].score == 0 && !first_pass){
                            if (j > previous_row_last_changed){
                                break innerloop;
                            }
                            else continue
                        }
                        matrix[i][j].score = 0
                        matrix[i][j].direction = dir.NONE
                        this_row_last_changed = j
                        break;
                    
                    case up_score:
                        if (matrix[i][j].score == up_score  && !first_pass){
                            if (j > previous_row_last_changed){
                                break innerloop;
                            }
                            else continue
                        }
                        
                        matrix[i][j].score = up_score
                        matrix[i][j].direction = dir.UP
                        this_row_last_changed = j
                        break;
                    
                    case left_score:
                        if (matrix[i][j].score == left_score  && !first_pass){
                            if (j > previous_row_last_changed){
                                break innerloop;
                            }
                            else continue
                        }

                        matrix[i][j].score = left_score
                        matrix[i][j].direction = dir.LEFT
                        this_row_last_changed = j
                        break;
                    
                    case diagonal_score:
                        if (matrix[i][j].score == diagonal_score && !first_pass){
                            if (j > previous_row_last_changed){
                                break innerloop;
                            }
                            else continue
                        }
                        matrix[i][j].score = diagonal_score
                        matrix[i][j].direction = dir.DIAGONAL
                        this_row_last_changed = j
                        break;
                }
                if (matrix[i][j].score >= threshold){
                    top_scores[[i,j].toString()] = matrix[i][j].score   //add coordinate pair to top_scores
                }
                else 
                    delete top_scores[[i,j].toString()]     //deletes coordinate from top_scores if it was there before
            }   
            if (first_pass)
                previous_row_last_changed = -1 //dont ever short circuit this calculation
            else
                previous_row_last_changed = this_row_last_changed;
            
        }
        
        var coords = Object.keys(top_scores)
        if (coords.length == 0) //no alignments pass the threshold
            return {}
        
        var max_score = top_scores[coords[0]]
        var max_coords = coords[0].split(',')
        for (var i=1; i < coords.length; i++) {
                if (top_scores[coords[i]] > max_score) {
                    max_score  = top_scores[coords[i]]
                    max_coords = coords[i].split(',')
                }
        }

        first_pass = false
        var n = parseInt(max_coords[0])
        var m = parseInt(max_coords[1])
        delete top_scores[[n,m].toString()]
        var seqv_result = ""
        var seqh_result = ""
        var path = []
        used_set[[n,m].toString()] = true

        while (matrix[n][m].direction != dir.NONE){
            path.push({row:n, col:m})
            switch( matrix[n][m].direction ){
                case dir.UP:
                    seqv_result = seqv[n] + seqv_result
                    seqh_result = "-" + seqh_result
                    n -= 1
                    break
                case dir.LEFT:
                    seqh_result = seqh[m] + seqh_result
                    seqv_result = "-" + seqv_result
                    m -= 1
                    break
                case dir.DIAGONAL:
                    used_set[[n,m].toString()] = true
                    seqv_result = seqv[n] + seqv_result
                    seqh_result = seqh[m] + seqh_result
                    n -= 1
                    m -= 1
                    break
            }
        }
        //Since gap penalties are negative, we can say for sure that the start point will be diagonal, so we start next calculation at n+1,m+1
        return {path: path, score: max_score, seqh_result: seqh_result, seqv_result: seqv_result, start_row:n+1, start_col:m+1}

    }
    
    var alignments = []
    var result = smith_waterman(1,1);
    if (Object.keys(result).length !=0)
        alignments.push(result)
    
    while (Object.keys(top_scores).length > 0) {
        result = smith_waterman(result.start_row,result.start_col);
        if (Object.keys(result).length !=0)
            alignments.push(result);

    }
       

    return alignments
}


function pprint_scores(matrix){
    return JSON.stringify(matrix.map(function(x){ return JSON.stringify(x.map(function(y){ return y.score}))})).split('","').join("\n").replace(/,/g," ").replace(/"|\[|\]/g,"")
}

	// body...

function k_best_scores(gapopen_pen, gapextend_pen, seqh, seqv, scoring_matrix,k) 
{

    
    //setting up scoring matrix and amino acid lookup
    var amino_acid_dict = {}
    var amino_acids = "A   R   N   D   C   Q  E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X".split(/\W+/)
    var i,j;
    for (i=1; i < amino_acids.length; i++) {
        amino_acid_dict[amino_acids[i-1]] = i
    }
    var score_matrix = get_scoring_matrix(scoring_matrix)

    //directions enum
    var dir = Object.freeze({NONE: 0, UP:1, LEFT:2, DIAGONAL:3}) // alignment direction
    
    //Prepend character to make string indexes line up with matrix indexes
    seqv = " " + seqv
    seqh = " " + seqh 
    if (gapopen_pen > 0)
        gapopen_pen = -gapopen_pen
    if (gapextend_pen > 0)
        gapextend_pen = -gapextend_pen
    //initialize matrix
    var matrix = new Array(seqv.length);
    
    for (i=0; i < matrix.length; i++) {
        matrix[i] = new Array(seqh.length);
    }

    for (i=0; i < matrix.length; i++) {
        for (j=0; j < matrix[0].length; j++) {
            matrix[i][j] = { direction: dir.NONE, score: 0 }
        }
    }
    var used_set = {}

    //the start point represents where you have to recalculate from i.e. the start of the last solution
    function smith_waterman(srow,scol) {  
        var up, left, diagonal, up_score, left_score, diagonal_score, max_score, max_coords
        var i
        max_score = 0
        max_coords = {row:0,col:0}
        
        
        
        for (i=1; i < matrix.length; i++) {
            for (j=1; j < matrix[0].length; j++) {
                if (i < srow || j < scol) {
                    //all we need from the cells to the top or left of last alignment are the max_score
                    if (matrix[i][j].score > max_score){
                        max_score = matrix[i][j].score
                        max_coords.row = i;
                        max_coords.col = j;    
                    }
                    continue;
                }

                up = matrix[i-1][j]
                left = matrix[i][j-1]
                diagonal = matrix[i-1][j-1]
                
                if (left.direction == dir.LEFT)
                    left_score = left.score + gapextend_pen
                else 
                    left_score = left.score + gapopen_pen
 
                if (up.direction == dir.UP)
                    up_score = up.score + gapextend_pen
                else
                    up_score = up.score + gapopen_pen
 
                if ([i,j].toString() in used_set)
                    diagonal_score = -999           // in Smith waterman
                else
                    diagonal_score = diagonal.score + score_matrix[amino_acid_dict[seqv[i]]][amino_acid_dict[seqh[j]]]
 
                switch (Math.max(diagonal_score,0,up_score,left_score)){
                    case 0:
                        matrix[i][j].score = 0
                        matrix[i][j].direction = dir.NONE
                        break;
                    case up_score:
                        matrix[i][j].score = up_score
                        matrix[i][j].direction = dir.UP
                        break;
                    case left_score:
                        matrix[i][j].score = left_score
                        matrix[i][j].direction = dir.LEFT
                        break;
                    case diagonal_score:
                        matrix[i][j].score = diagonal_score
                        matrix[i][j].direction = dir.DIAGONAL
                        break;
                }
                if (matrix[i][j].score > max_score){
                    max_score = matrix[i][j].score
                        max_coords.row = i;
                        max_coords.col = j;    
                }

            }   
            
        }
        

        //Collect the maximum score. Unfortunately we go over the entire matrix for now
        
    
        var n = max_coords.row
        var m = max_coords.col
        var seqv_result = ""
        var seqh_result = ""
        var path = []
        used_set[[n,m].toString()] = true
        while (matrix[n][m].direction != dir.NONE){
            path.push({row:n, col:m})
            switch( matrix[n][m].direction ){
                case dir.UP:
                    seqv_result = seqv[n] + seqv_result
                    seqh_result = "-" + seqh_result
                    n -= 1
                    break
                case dir.LEFT:
                    seqh_result = seqh[m] + seqh_result
                    seqv_result = "-" + seqv_result
                    m -= 1
                    break
                case dir.DIAGONAL:
                    used_set[[n,m].toString()] = true
                    seqv_result = seqv[n] + seqv_result
                    seqh_result = seqh[m] + seqh_result
                    n -= 1
                    m -= 1
                    break
            }
        }
        //If gap penalties are always negative, we can say for sure that the end point will be diagonal, so we start next calculation at n+1,m+1
        return {path: path, score: max_score, seqh_result: seqh_result, seqv_result: seqv_result, start_row:n+1, start_col:m+1}

    }
    var result = smith_waterman(1,1)
    
    
    var alignments = [result]
    for (i=1; i < k; i++) {
        result = smith_waterman(result.start_row,result.start_col)
        alignments.push(result)

    }
       

    return alignments
}


function pprint_scores(matrix){
    return JSON.stringify(matrix.map(function(x){ return JSON.stringify(x.map(function(y){ return y.score}))})).split('","').join("\n").replace(/,/g," ").replace(/"|\[|\]/g,"")
}

    // body...

