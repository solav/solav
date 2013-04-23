

tests = {};
//freuencies From http://web.expasy.org/protscale/pscale/A.A.Swiss-Prot.html Bairoch A.
tests.frequencies =[0,8.25, 5.5, 4.06, 5.45, 1.37, 3.93, 6.75, 7.07, 2.27, 5.96, 9.66, 5.84, 2.42, 3.86, 4.70, 6.5, 5.34, 1.08, 2.92, 6.87] ;
tests.total = 0;
tests.cumulative_frequencies = tests.frequencies.map(function(x){return tests.total+=x;});


var amino_acids = "A   R   N   D   C   Q  E   G   H   I   L   K   M   F   P   S   T   W   Y   V  ".split(/\W+/);


function randomString(n){
	var string = '';
	for (var i =0; i < n; i++ ){
        var val = Math.random()*tests.total;
        for (var j =0; tests.cumulative_frequencies[j]<val;j++){      	
        }
     	string += amino_acids[j-1];
	}
    return string;
}

function sequencelength_speedTest(size){
    //unrecorded one to get the optimization going, otherwise the first ones take a bizarrely long time
    var s1 = randomString(10)
    var s2 = randomString(10)
    var results = k_best_scores(10,5,s1,s2,"BLOSUM62", 1)

	for (var i=0; i<10; i++){		
        var s1 = randomString(size)
        var s2 = randomString(size)

        console.time("k_best " + size)
        var results = k_best_scores(10,5,s1,s2,"BLOSUM62", Math.min(size,10))
        console.timeEnd("k_best " + size)
        
        threshhold = results[results.length-1].score

        console.time("watermaneggert " + size)
        var results = waterman_eggert(10,5,s1,s2,"BLOSUM62",threshhold)
        console.timeEnd("watermaneggert " + size)
    }
	
}

function k_speedTest(k) {

    //how does changing k effect each algorithm

    var s1 = randomString(10)
    var s2 = randomString(10)
    var results = k_best_scores(10,5,s1,s2,"BLOSUM62", 1)

    for (var i=0; i<10; i++){       
        var s1 = randomString(500)
        var s2 = randomString(500)
        console.time("k_best " + size)
        var results = k_best_scores(10,5,s1,s2,"BLOSUM62", k)
        console.timeEnd("k_best " + size)
        
        threshhold = results[results.length-1].score

        console.time("watermaneggert " + size)
        var results = waterman_eggert(10,5,s1,s2,"BLOSUM62",threshhold)
        console.timeEnd("watermaneggert " + size)
    }
    
}    

	

function speedTest(){
	sizes = [10,20,50,100,200,300,400,500,1000,2000]
	//sizes = [10000]
	for (var j=0; j<sizes.length; j++)
		speedTest(sizes[j])
}


// obsolete first version 
// function watermaneggert_basic(gapopen_pen, gapextend_pen, seqh, seqv, scoring_matrix,k) {  
//     //setting up scoring matrix and amino acid lookup
//     var amino_acid_dict = {}
//     var amino_acids = "A   R   N   D   C   Q  E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X".split(/\W+/)
//     var i,j;
//     for (i=1; i < amino_acids.length; i++) {
//         amino_acid_dict[amino_acids[i-1]] = i
//     }
//     var score_matrix = get_scoring_matrix(scoring_matrix)

//     //directions enum
//     var dir = Object.freeze({NONE: 0, UP:1, LEFT:2, DIAGONAL:3}) // alignment direction
    
//     //Prepend character to make string indexes line up with matrix indexes
//     seqv = " " + seqv
//     seqh = " " + seqh 
    
//     //initialize matrix
//     var matrix = new Array(seqv.length);
    
//     for (i=0; i < matrix.length; i++) {
//         matrix[i] = new Array(seqh.length);
//     }

//     for (i=0; i < matrix.length; i++) {
//         for (j=0; j < matrix[0].length; j++) {
//             matrix[i][j] = { direction: dir.NONE,  score: 0 }
//         }
//     }
//     var used_set = {}

//     //the start point represents where you have to recalculate from i.e. the start of the last solution

//     function smith_waterman(srow,scol) {  
//         var up, left, diagonal, up_score, left_score, diagonal_score, max_score, max_coords
//         var iZ
        
        
//         for (i=srow; i < matrix.length; i++) {
//             for (j=scol; j < matrix[0].length; j++) {
//                 up = matrix[i-1][j]
//                 left = matrix[i][j-1]
//                 diagonal = matrix[i-1][j-1]
                
//                 if (left.left_gap)
//                     left_score = left.score - gapextend_pen
//                 else 
//                     left_score = left.score - gapopen_pen
 
//                 if (up.up_gap)
//                     up_score = up.score - gapextend_pen
//                 else
//                     up_score = up.score - gapopen_pen
 
//                 if ([i,j].toString() in used_set)
//                     diagonal_score = -999           // in Smith waterman
//                 else
//                     diagonal_score = diagonal.score + score_matrix[amino_acid_dict[seqv[i]]][amino_acid_dict[seqh[j]]]
 
//                 switch (Math.max(diagonal_score,0,up_score,left_score)){
//                     case 0:
//                         matrix[i][j].score = 0
//                         matrix[i][j].up_gap = false
//                         matrix[i][j].left_gap = false
//                         matrix[i][j].direction = dir.NONE
//                         break;
//                     case up_score:
//                         matrix[i][j].score = up_score
//                         matrix[i][j].up_gap = true
//                         matrix[i][j].left_gap = false
//                         matrix[i][j].direction = dir.UP
//                         break;
//                     case left_score:
//                         matrix[i][j].score = left_score
//                         matrix[i][j].up_gap = false
//                         matrix[i][j].left_gap = true
//                         matrix[i][j].direction = dir.LEFT
//                         break;
//                     case diagonal_score:
//                         matrix[i][j].score = diagonal_score
//                         matrix[i][j].up_gap = false
//                         matrix[i][j].left_gap = false
//                         matrix[i][j].direction = dir.DIAGONAL
//                         break;
//                 }
//             }   
            
//         }
        

//         //Collect the maximum score. Unfortunately we go over the entire matrix for now
//         max_score = 0
//         max_coords = {row:0,col:0}
//         for (var i=1; i < matrix.length; i++) {
//             for (var j=1; j < matrix[0].length; j++) {
//                 if (matrix[i][j].score > max_score){
//                     max_score = matrix[i][j].score
//                     max_coords.row = i;
//                     max_coords.col = j;
//                 }
//             }
//         }
        
    
//         var n = max_coords.row
//         var m = max_coords.col
//         var seqv_result = ""
//         var seqh_result = ""
//         var path = []
//         used_set[[n,m].toString()] = true
//         while (matrix[n][m].score != dir.NONE){
            
//             path.push({row:n, col:m})
//             switch( matrix[n][m].direction ){
//                 case dir.UP:
//                     seqv_result = seqv[n] + seqv_result
//                     seqh_result = "-" + seqh_result
//                     n -= 1
//                     break
//                 case dir.LEFT:
//                     seqh_result = seqh[m] + seqh_result
//                     seqv_result = "-" + seqv_result
//                     m -= 1
//                     break
//                 case dir.DIAGONAL:
//                     used_set[[n,m].toString()] = true
//                     seqv_result = seqv[n] + seqv_result
//                     seqh_result = seqh[m] + seqh_result
//                     n -= 1
//                     m -= 1
//                     break
				
//             }


//         }
//         //If gap penalties are always negative, we can say for sure that the end point will be diagonal, so we start next calculation at n+1,m+1
//         return {path: path, score: max_score, seqh_result: seqh_result, seqv_result: seqv_result, start_row:n+1, start_col:m+1}

//     }
//     var result = smith_waterman(1,1)
    
//     var alignments = [result]
//     for (i=1; i < k; i++) {
//         console.log(result.a1)
//         console.log(result.a2)
//         result = smith_waterman(result.start_row,result.start_col)
//         alignments.push(result)

//     }
//     console.log(result.a1)
//     console.log(result.a2)
//     return alignments
// }
// //for debug
