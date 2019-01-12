(*********************************************************************************************************)
(****************************************** initialization ***********************************************)
(*********************************************************************************************************)

(************************************* boundary type: Random[0, d] ***************************************)
(*
This function initializes the boundary values of the given matrix to generate the desire initial surface.

Input:
		inMatrix_: a zero matrix.
    randomRange_: the cap of the RandomReal for the outer boundary values.
    innerBoundaryRandomRange: the cap of the RandomReal for the inner boundary\[NonBreakingSpace]values.
	  innerBoundary_: boolean variable for inner boundary.
	  innerBoundaryXIndex_: the x index of the top right corner of the inner boundary.
	  innerBoundaryYIndex_: the y index of the top right corner of the inner boundary.
	  innerBoundaryXwidth_: the width of the inner boundary.
	  innerBoundaryYwidth_: the height of the inner boundary.
	  randomRangeDomain_: the cap of the RandomReal for the interior lattice points.

Output:
    The initialized matrix
*)
randomInitilization[inMatrix_, randomRange_, innerBoundaryRandomRange_, innerBoundary_, innerBoundaryXIndex_, 
                    innerBoundaryYIndex_, innerBoundaryXwidth_, innerBoundaryYwidth_, randomRangeDomain_] := 
   	Block[
        {m = inMatrix, h = Dimensions[inMatrix][[1]], l = Dimensions[inMatrix][[2]], r = innerBoundary, 
        rx = innerBoundaryXIndex, ry = innerBoundaryYIndex, rxl = innerBoundaryXwidth, 
        ryl = innerBoundaryYwidth, rv = randomRange, rrv = innerBoundaryRandomRange, 
        rd = randomRangeDomain},

      	m = Table[0, {x, 1 , h}, {y, 1, l}];
      	If[rx >= h, rx = Floor[h / 3]];
      	If[ry >= l, ry = Floor[l/ 3]];
      	If[rx + rxl >= h, rxl = h - 1 - rx];
      	If[ry + ryl >= l, ryl = l - 1 - ry];
      	If[r,
      			Do[
                If[x == 1 || x == h || y == 1 || y == l, 
                    m[[x, y]] = RandomReal[rv], m[[x, y]] = RandomReal[rd]
                ],
                {y, 1, l}, {x, 1, h}
            ];
      			Do[
                If[x == rx || x == rx + rxl || y == ry || y == ry + ryl, 
                    m[[x, y]] = RandomReal[rrv]
                ],
                {x, rx, rx + rxl}, {y, ry, ry + ryl}
            ];
      			,
      			Do[
                If[x == 1 || x == h || y == 1 || y == l, 
                    m[[x, y]] = RandomReal[rv], m[[x, y]] = RandomReal[rd]
                ],
                {y, 1, l}, {x, 1, h}
            ];
        ];
        m
    ]

(************************************* boundary type: c + sin(ax + by) ***************************************)
(*
Input:           
		inMatrix_: a zero matrix.
	  innerBoundary_: boolean variable for inner boundary.
	  innerBoundaryXIndex_: the x index of the top right corner of the inner boundary
	  innerBoundaryYIndex_: the y index of the top right corner of the inner boundary
	  innerBoundaryXwidth_: the width of the inner boundary.
	  innerBoundaryYwidth_: the height of the inner boundary.
	  inA_: the coefficient a of the outer boundary.
	  inB_: the coefficient b of the outer boundary.
	  inC_: the coefficient c of the outer boundary.
	  innerBoundaryA_: the coefficient a of the inner boundary.
	  innerBoundaryB_: the coefficient b of the inner boundary.
	  innerBoundaryC_: the coefficient c of the inner boundary.
	  randomRangeSinDomain_: the cap of the RandomReal for the interior lattice points.
Output:
    The initialized matrix.
*)
sinInitialization[inMatrix_, innerBoundary_, innerBoundaryXIndex_, innerBoundaryYIndex_, 
                  innerBoundaryXwidth_, innerBoundaryYwidth_, inA_, inB_, inC_, 
                  innerBoundaryA_, innerBoundaryB_, innerBoundaryC_, randomRangeSinDomain_] := 
    Block[
        {m = inMatrix, h = Dimensions[inMatrix][[1]], l = Dimensions[inMatrix][[2]], 
        r = innerBoundary, rx = innerBoundaryXIndex, ry = innerBoundaryYIndex, 
        rxl = innerBoundaryXwidth, ryl = innerBoundaryYwidth, a = inA, b = inB, 
        c = inC, ra = innerBoundaryA, rb = innerBoundaryB, rc = innerBoundaryC, 
        rds = randomRangeSinDomain},

      	m = Table[0, {x, 1 , h}, {y, 1, l}];
      	If[rx >= h, rx = Floor[h / 3]];
      	If[ry >= l, ry = Floor[l/ 3]];
      	If[rx + rxl >= h, rxl = h - 1 - rx];
      	If[ry + ryl >= l, ryl = l - 1 - ry];
      	If[r,
      			Do[
                If[x == 1 || x == h || y == 1 || y == l, 
                    m[[x, y]] = Sin[a*x + b*y] + c, 
                    m[[x, y]] = RandomReal[{-rds, rds}]
                ],
                {y, 1, l}, {x, 1, h}
            ];
      			Do[
                If[x == rx || x == rx + rxl || y == ry || y == ry + ryl, 
                    m[[x, y]] = Sin[ra*x + rb*y] + rc
                ],
                {x, rx, rx + rxl}, {y, ry, ry + ryl}
            ];
      			, 
      			Do[
                If[x == 1 || x == h || y == 1 || y == l, 
                    m[[x, y]] = Sin[a*x + b*y] + c, 
                    m[[x, y]] = RandomReal[{-rds, rds}]
                ],
                {y, 1, l}, {x, 1, h}
            ];
        ];
        m
    ]

(************************************* boundary type: c + a x^2 + b y^2 ***************************************)
(*
Input:
		inMatrix_: a zero matrix.
	  innerBoundary_: boolean variable for inner boundary.
	  innerBoundaryXIndex_: the x index of the top right corner of the inner boundary.
	  innerBoundaryYIndex_: the y index of the top right corner of the inner boundary.
	  innerBoundaryXwidth_: the width of the inner boundary.
	  innerBoundaryYwidth_: the height of the inner boundary.
	  inA_: the coefficient a of the outer boundary.
	  inB_: the coefficient b of the outer boundary.
	  inC_: the coefficient c of the outer boundary.
	  innerBoundaryA_: the coefficient a of the inner boundary.
	  innerBoundaryB_: the coefficient b of the inner boundary.
    innerBoundaryC_: the coefficient c of the inner boundary.
	  randomRangeQuadraDomain_: the cap of the RandomReal for the interior lattice points.
Output:
    The initialized matrix.
*)
quadraticInitialization[inMatrix_, innerBoundary_, innerBoundaryXIndex_, innerBoundaryYIndex_, 
                        innerBoundaryXwidth_, innerBoundaryYwidth_, inA_, inB_, inC_, 
                        innerBoundaryA_, innerBoundaryB_ , innerBoundaryC_, randomRangeQuadraDomain_] :=
    Block[
        {m = N[inMatrix], h = Dimensions[inMatrix][[1]], l = Dimensions[inMatrix][[2]], 
        r = innerBoundary, rx = innerBoundaryXIndex, ry = innerBoundaryYIndex, 
        rxl = innerBoundaryXwidth, ryl = innerBoundaryYwidth, a = inA, b = inB, c = inC, 
        ra = innerBoundaryA, rb = innerBoundaryB, rc = innerBoundaryC, rdq = randomRangeQuadraDomain},
    	
        m = Table[0, {x, 1 , h}, {y, 1, l}];
      	If[rx >= h, rx = Floor[h / 3]];
      	If[ry >= l, ry = Floor[l/ 3]];
      	If[rx + rxl >= h, rxl = h - 1 - rx];
      	If[ry + ryl >= l, ryl = l - 1 - ry];
      	If[r,
      			Do[
                If[x == 1 || x == h || y == 1 || y == l, 
                    m[[x, y]] = N[a*(x^2)  + b*(y^2) + rc * x * y], 
                    m[[x, y]] = RandomReal[{-rdq, rdq}]
                ],
                {y, 1, l}, {x, 1, h}
            ];
      			Do[
                If[x == rx || x == rx + rxl || y == ry || y == ry + ryl, 
                    m[[x, y]] = ra*(x^2)  + rb*(y^2) + rc * x * y
                ],
                {x, rx, rx + rxl}, {y, ry, ry + ryl}
            ];
      			, 
      			Do[
                If[x == 1 || x == h || y == 1 || y == l, 
                    m[[x, y]] = N[a*(x^2)  + b*(y^2) + rc * x * y], 
                    m[[x, y]] = RandomReal[{-rdq, rdq}]
                ],
                {y, 1, l}, {x, 1, h}
            ];
        ];
        m
    ]

(*********************************************************************************************************)
(**************************************** relaxation process *********************************************)
(*********************************************************************************************************)

(*
This function applies one step of the relaxation process to a given matrix and output the resulting matrix.

Input:
    inMatrix_: he matrix we want to apply method of relaxation on.
Output:
    Result matrix after one iteration.
*)
relaxationProcess[inMatrix_] := 
    Block[
        {x = N[inMatrix], h = Dimensions[inMatrix][[1]], l = Dimensions[inMatrix][[2]]},

        Do[
            x[[i, j]] = (x[[i + 1, j]] + x[[i - 1, j]] + x[[i, j + 1]] + x[[i, j - 1]])/4, 
            {j , 2 , l - 1}, {i, 2, h - 1}
        ];
        x[[1, 1]] = (x[[2, 1]] + x[[1, 2]]) / 2;
        x[[1, l]] = (x[[2, l]] + x[[1, l - 1]]) / 2;
        x[[h, 1]] = (x[[h - 1, 1]] + x[[h, 2]]) / 2;
        x[[h, l]] = (x[[h - 1, l]] + x[[h, l - 1]]) / 2;
        x
    ]

(************************************* case with inner boundary ***************************************)
(*
Input:
    inMatrix_: the matrix we want to apply method of relaxation on.
    ininnerBoundaryX_: the x index of the top right corner of the inner boundary.
    ininnerBoundaryY_: the y index of the top right corner of the inner boundary.
    ininnerBoundaryXLen_: the width of the inner boundary.
    ininnerBoundaryYLen_: the height of the inner boundary.

Output:
    Result matrix after one iteration.
*)
relaxationProcessInnerBoundary[inMatrix_, ininnerBoundaryX_, ininnerBoundaryY_, 
                               ininnerBoundaryXLen_, ininnerBoundaryYLen_] := 
    Block[
        {m = N[inMatrix], h = Dimensions[inMatrix][[1]], l = Dimensions[inMatrix][[2]], 
        x = ininnerBoundaryX, y = ininnerBoundaryY, xl = ininnerBoundaryXLen, 
        yl = ininnerBoundaryYLen},

        If[x >= h, x = Floor[h / 3]];
        If[y >= l, y = Floor[l/ 3]];
        If[x + xl >= h, xl = h - 1 - x];
        If[y + yl >= l, yl = l - 1 - y];
        Do[
            m[[i, j]] = (m[[i + 1, j]] + m[[i - 1, j]] + m[[i, j + 1]] + 
            m[[i, j - 1]])/4, 
            {j , 2 , l - 1}, {i, 2, x - 1}
        ]; (*above inner boundary*)
        Do[
            m[[i, j]] = (m[[i + 1, j]] + m[[i - 1, j]] + m[[i, j + 1]] + 
            m[[i, j - 1]])/4, 
            {j , 2 , l - 1}, {i, x + xl + 1, h - 1}
        ]; (*below inner boundary*)
        Do[
            m[[i, j]] = (m[[i + 1, j]] + m[[i - 1, j]] + m[[i, j + 1]] + 
            m[[i, j - 1]])/4, 
            {j , 2 , y - 1}, {i, x, x + xl}
        ]; (*left of inner boundary*)
        Do[
            m[[i, j]] = (m[[i + 1, j]] + m[[i - 1, j]] + m[[i, j + 1]] + 
            m[[i, j - 1]])/4, {j , y + 1 + yl, l - 1}, {i, x, x + xl}
        ]; (*right of inner boundary*)
        m
    ]

(*
This function assigns values to the four corners of the matrix as the mean value of the two sourrounding points.

Input: 
    inMatrix_: the Matrix whose four corner points are unassigned.
    currentIteration_: the number of iteration.
Output: 
    The matrix with all latice points assigned with correct values.
*)
initializeCorner[inMatrix_, currentIteration_] := 
    Block[
        {m = N[inMatrix], h = Dimensions[inMatrix[[1]]][[1]], 
        l = Dimensions[inMatrix[[1]]][[2]], c = currentIteration},

        m[[c, 1, 1]] = (m[[c, 1, 2]] + m[[c, 2, 1]]) / 2;
        m[[c, 1, l]] = (m[[c, 1, l - 1]] + m[[c, 2, l]]) / 2;
        m[[c, h, 1]] = (m[[c, h - 1, 1]] + m[[c, h, 2]]) / 2;
        m[[c, h, l]] = (m[[c, h - 1, l]] + m[[c, h, l - 1]]) / 2;
        m
    ]

(*********************************************************************************************************)
(***************************************** Dirichlet energy **********************************************)
(*********************************************************************************************************)

(*
This function calculates the Dirichlet energy of given matrix
	
Input: 
    inMatrix_: the matrix representing the surface we want to calculate the energy on.
Output: 
    The Dirichlet energy of the given matrix.
*)
dirichletEnergy[inMatrix_] := 
    Block[
        {m = N[inMatrix], h = Dimensions[inMatrix][[1]], l = Dimensions[inMatrix][[2]], 
        sum = 0},

        Do[
            sum += (m[[x, y]] - m[[x, y - 1]])^2 + (m[[x, y]] - m[[x - 1, y]])^2, 
            {y, 2, l}, {x, 2, h}
        ];
        Do[
            sum += (m[[x , 1]] - m[[x + 1, 1]])^2, 
            {x, 1, h - 1}
        ];
        Do[
            sum += (m[[1, y]] - m[[1, y + 1]])^2, {y, 1, l - 1}
        ];
        sum
    ]

(************************************* case with inner boundary ***************************************)
(*	
Input:
    inMatrix_: the matrix representing the surface we want to calculate the energy on.
    innerBoundaryXIndex_: the x index of the top right corner of the inner boundary.
    innerBoundaryYIndex_: the y index of the top right corner of the inner boundary.
    innerBoundaryXwidth_: the width of the inner boundary.
    innerBoundaryXwidth_: the height of the inner boundary.
Output: 
    The Dirichlet energy of the given matrix.
*)
dirichletEnergyInnerBoundary[inMatrix_, innerBoundaryXIndex_, innerBoundaryYIndex_, 
                             innerBoundaryXwidth_, innerBoundaryYwidth_] := 
    Block[
        {m = N[inMatrix], h = Dimensions[inMatrix][[1]], l = Dimensions[inMatrix][[2]], 
        rX = innerBoundaryXIndex, rXL = innerBoundaryXwidth, rY = innerBoundaryYIndex, 
        rYL = innerBoundaryYwidth, sum = 0},

        If[rX >= h, rX = Floor[h / 3]];
        If[rY >= l, rY = Floor[l/ 3]];
        If[rX + rXL >= h, rXL = h - 1 - rX];
        If[rY + rYL >= l, rYL = l - 1 - rY];
        Do[
            sum += (m[[x, y]] - m[[x, y - 1]])^2 + (m[[x, y]] - m[[x - 1, y]])^2, 
            {y, 2, l}, {x, 2, rX }
        ];
        Do[
            sum += (m[[x, y]] - m[[x, y - 1]])^2 + (m[[x, y]] - m[[x - 1, y]])^2,
            {y, 2, l}, {x, rX + rXL , h}
        ];
        Do[
            sum += (m[[x, y]] - m[[x, y - 1]])^2 + (m[[x, y]] - m[[x - 1, y]])^2, 
            {y, 2, rY}, {x, rX, rX + rXL}
        ];
        Do[
            sum += (m[[x, y]] - m[[x, y - 1]])^2 + (m[[x, y]] - m[[x - 1, y]])^2, 
            {y, rY + rYL, l}, {x, rX, rX + rXL}
        ];
        Do[
            sum += (m[[x , 1]] - m[[x + 1, 1]])^2, 
            {x, 1, h - 1}
        ];
        Do[
            sum += (m[[1, y]] - m[[1, y + 1]])^2, 
            {y, 1, l - 1}
        ];
        Do[
            sum += (m[[x , 1]] - m[[x + 1, 1]])^2, 
            {x, rX, rX + rXL - 1}
        ];
        Do[
            sum += (m[[1, y]] - m[[1, y + 1]])^2, 
            {y, rY, rY + rYL - 1}
        ];
        sum
    ]

(*********************************************************************************************************)
(**************************************** ploting function **********************************************)
(*********************************************************************************************************)

(*
This function is used to plot the given surface and also displays the energy of the surface.

Input:    
    inMatrix_: the matrix representing the surface we want to calculate the energy on.
    energy_: the value of the Dirichlet energy we want to label on the graph.
    boundaryType_: the type of boundary function we want to plot(random, trigonometric, or quadratic).
    randomRange_: the number of random value range we want to display on the graph.
    sinA_: the coefficient a for the trigonometric boundary function.
    sinB_: the coefficient b for the trigonometric boundary function.
    sinC_: the coefficient c for the trigonometric boundary function.
    quaA_: the coefficient a for the quadratic boundary function.
    quaB_: the coefficient b for the quadratic boundary function.
    quaC_: the coefficient c for the quadratic boundary function. 
    currentIteration_: the number of iteration.

Output:
    The graph of the discrete harmonic function surface.	
*)     
plotGraph[inMatrix_, energy_, boundaryType_, randomRange_, sinA_, sinB_, 
          sinC_, quaA_, quaB_, quaC_, currentIteration_] := 
    Block[
        {m = N[inMatrix], f = energy, h = Dimensions[inMatrix][[1]], 
        l = Dimensions[inMatrix][[2]], g = boundaryType, r = randomRange, 
        qa = quaA, qb = quaB, qc = quaC, cs = currentIteration, sa = sinA, 
        sb = sinB, sc = sinC}, 
        Switch[g,
            "Random[0, d]",	
                Column[{
                    Column[{
                        Row[{
                            Style["Outer boundary: Random[ 0, ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                            Style[r, FontSize -> 12] Style["]", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"]
                        }], 
                        Row[{
                            Style["Number of iterations: ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                            Style[cs, FontSize -> 12]
                        }]
                    }, Alignment -> Center], 
                    ListPlot3D[m, Mesh -> {Min[l, 15], Min[h, 15]},  ColorFunction -> "CMYKColors", 
                                  PlotStyle -> Opacity[0.5], ImageSize -> {300, 300}, 
                                  SphericalRegion -> True, Filling -> Bottom, PlotRange -> All, 
                                  ImagePadding -> 15, 
                                  PlotLabel -> 
                                      Row[{
                                          Style["Dirichlet Energy: ", Italic,  FontSize -> 12, 
                                          FontFamily -> "Arial Hebrew Scholar"], 
                                          Style[f, FontSize -> 12]
                                      }], 
                    PerformanceGoal -> "Speed"]
                }, Alignment -> Center],	
            "c + sin(ax + by)",
                Column[{
                    Column[{
                        Row[{
                            Style["Outer boundary: ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                            Style[StringJoin[ToString[sc], " + ", "sin(", ToString[sa], "x + ", 
                                  ToString[sb], "y)"], FontSize -> 12]
                        }],
                        Row[{
                            Style["Number of iterations: ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                            Style[cs, FontSize -> 12]
                        }]
                    }, Alignment -> Center], 
                    ListPlot3D[m, Mesh -> {Min[l, 15], Min[h, 15]},  ColorFunction -> "BrightBands", 
                                  PlotStyle -> Opacity[0.5], ImageSize -> {300, 300}, 
                                  SphericalRegion -> True, Filling -> Bottom, PlotRange -> All, 
                                  ImagePadding -> 15, 
                                  PlotLabel -> 
                                      Row[{
                                          Style["Dirichlet Energy: ", Italic,  FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                                          Style[f, FontSize -> 12]
                                      }], 
                    PerformanceGoal -> "Speed"]
                }, Alignment -> Center],
            TraditionalForm[c + (ax^2 + bx^2)],
                Column[{
                    Column[{
                        Row[{
                            Style["Outer boundary: ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                            Style[ StringJoin[ ToString[qa], "x^2 + ", ToString[qb], "y^2 "], FontSize -> 12]
                        }],
                    Row[{
                        Style["Number of iterations: ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                        Style[cs, FontSize -> 12]
                    }]
                }, Alignment -> Center], 
                ListPlot3D[m, Mesh -> {Min[l, 15], Min[h, 15]},  ColorFunction -> "LightTemperatureMap",
                              PlotStyle -> Opacity[0.5], ImageSize -> {300, 300}, SphericalRegion -> True, 
                              Filling -> Bottom, PlotRange -> All, ImagePadding -> 15, 
                              PlotLabel -> 
                                  Row[{
                                      Style["Dirichlet Energy: ", Italic,  FontSize -> 12, FontsFamily -> "Arial Hebrew Scholar"], 
                                      Style[f,  FontSize -> 12]
                                  }], 
                PerformanceGoal -> "Speed"]
            }, Alignment -> Center]
        ]
    ]

(************************************* case with inner boundary ***************************************)
(*
Input: 
    inMatrix_: the matrix representing the surface we want to calculate the energy on.
    energy_: the value of the Dirichlet energy we want to label on the graph.
    innerBoundaryXIndex_: the x index of the top right corner of the inner boundary.
    innerBoundaryXwidth_: the width of the inner boundary.
    innerBoundaryYIndex_: the y index of the top right corner of the inner boundary.
    innerBoundaryYwidth_: the height of the inner boundary.
    boundaryType_: the type of boundary function we want to plot(random, trigonometric, or quadratic).
    randomRange_: the number of random value range we want to display on the graph.
    innerBoundaryRandomRange_: the number of random value range for the inner boundary we want to display on the graph.
    sinA_: the coefficient a for the trigonometric boundary function.
    sinB_: the coefficient b for the trigonometric boundary function.
    sinC_: the coefficient c for the trigonometric boundary function.
    innerBoundaryA_: the coefficient a for the inner trigonometric function.
    innerBoundaryB_: the coefficient b for the inner trigonometric function.
    innerBoundaryC_: the coefficient c for the inner trigonometric function.
    quaA_: the coefficient a for the quadratic boundary function.
    quaB_: the coefficient b for the quadratic boundary function.
    quaC_: the coefficient c for the quadratic boundary function.
    quaInnerBoundaryA_: the coefficient a for the inner quadratic function.
    quaInnerBoundaryB_: the coefficient b for the inner quadratic function.
    quaInnerBoundaryC_: the coefficient c for the inner quadratic function.
    currentIteration_: the number of iteration.

Output:
    The graph of the discrete harmonic function.
*)
plotGraphInnerBoundary[inMatrix_, energy_, innerBoundaryXIndex_, innerBoundaryXwidth_, 
                      innerBoundaryYIndex_, innerBoundaryYwidth_, boundaryType_, randomRange_, 
                      innerBoundaryRandomRange_, sinA_, sinB_, sinC_, innerBoundaryA_, 
                      innerBoundaryB_, innerBoundaryC_, quaA_, quaB_, quaC_, quaInnerBoundaryA_,
                      quaInnerBoundaryB_, quaInnerBoundaryC_, currentIteration_] :=

    Block[
        {m = inMatrix, l = Dimensions[inMatrix][[2]], h = Dimensions[inMatrix][[1]], 
        f = energy, rX = innerBoundaryXIndex, rXL = innerBoundaryXwidth, 
        rY = innerBoundaryYIndex, rYL = innerBoundaryYwidth, g = boundaryType,
        r = randomRange, rr = innerBoundaryRandomRange, ra = innerBoundaryA, 
        rb = innerBoundaryB, rc = innerBoundaryC, qa = quaA, qb = quaB, qc = quaC, 
        qra = quaInnerBoundaryA, qrb = quaInnerBoundaryB, qrc = quaInnerBoundaryC, 
        cs = currentIteration, sa = sinA, sb = sinB, sc = sinC}, 

        Switch[g,
            "Random[0, d]",
                Column[{
                    Column[{
                        Row[{
                            Style["Outer boundary: Random[ 0, ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                            Style[r, FontSize -> 12], Style["]", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"]
                        }],
                        Row[{
                            Style["Inner boundary: Random[ 0, ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                            Style[rr, FontSize -> 12], 
                            Style["]", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"]
                        }], 
                        Row[{
                            Style["Number of iterations: ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                            Style[cs, FontSize -> 12]
                        }]
                    }, Alignment -> Center], 
                    ListPlot3D[m, ColorFunction -> "CMYKColors", Mesh -> {Min[l, 15], Min[h, 15]}, 
                                  PlotStyle -> Opacity[0.5], Filling -> Bottom, PlotRange -> All, 
                                  ImageSize -> {300, 300}, SphericalRegion -> True, ImagePadding -> 15, 
                                  PlotLabel -> 
                                  Row[{
                                      Style["Dirichlet Energy: ", Italic,  FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                                      Style[f,  FontSize -> 12]}],
                                  RegionFunction -> 
                                      Function[{y, x}, 
                                            (1 <= x <= rX && 1 <= y <= l) (*inner boundary above*)
                                         || (rX + rXL <= x <= h && 1 <= y <= l) (*inner boundary below*)
                                         || ((rX <= x <= rX + rXL) && 1 <= y <= rY)  (*inner boundary left*)
                                         || ((rX <= x <= rX + rXL) && rY + rYL + 1 <= y <= l )  (*inner boundary right*)
                                      ], 
                              PerformanceGoal -> "Speed"]
                }, Alignment -> Center],
            "c + sin(ax + by)",
                Column[{
                    Column[{
                        Row[{
                            Style["Outer boundary: ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                            Style[StringJoin[ ToString[sc], " + ", "sin(", ToString[sa], 
                            "x + ", ToString[sb], "y)"], FontSize -> 12]
                        }], 
                        Row[{
                            Style["Inner boundary: ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                            Style[StringJoin[ToString[rc], " + ", "sin(", ToString[ra], "x + ", ToString[rb], "y)"], FontSize -> 12]
                        }], 
                        Row[{
                            Style["Number of iterations: ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                            Style[cs, FontSize -> 12]
                        }]
                    }, Alignment -> Center], 
                    ListPlot3D[m, ColorFunction -> "BrightBands", Mesh -> {Min[l, 15], Min[h, 15]}, 
                                  PlotStyle -> Opacity[0.5], Filling -> Bottom, PlotRange -> All, 
                                  ImageSize -> {300, 300}, SphericalRegion -> True, ImagePadding -> 15, 
                                  PlotLabel -> 
                                      Row[{
                                          Style["Dirichlet Energy: ", Italic,  FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                                          Style[f,  FontSize -> 12]
                                      }],
                                  RegionFunction -> 
                                      Function[{y, x}, 
                                          (1 <= x <= rX && 1 <= y <= l) (*inner boundary above*)
                                       || (rX + rXL <= x <= h && 1 <= y <= l) (*inner boundary below*)
                                       || ((rX <= x <= rX + rXL) && 1 <= y <= rY) (*inner boundary left*)
                                       || ((rX <= x <= rX + rXL) && rY + rYL + 1 <= y <= l) (*inner boundary right*)], 
                    PerformanceGoal -> "Speed"]
                }, Alignment -> Center],
            TraditionalForm[c + (ax^2 + bx^2)],
                Column[{
                    Column[{
                        Row[{
                            Style["Outer boundary: ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                            Style[ StringJoin[ToString[qc], "+ (", ToString[qa], "x^2 + ", ToString[qb], "y^2)"], FontSize -> 12]
                        }], 
                    Row[{
                        Style["Inner boundary: ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                        Style[StringJoin[ ToString[qrc], "+ (", ToString[qra], "x^2 + ", ToString[qrb], "y^2)"], FontSize -> 12]
                    }], 
                    Row[{
                        Style["Number of iterations: ", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                        Style[cs, FontSize -> 12]
                    }]
                }, Alignment -> Center], 
                ListPlot3D[m, ColorFunction -> "LightTemperatureMap", 
                              Mesh -> {Min[l, 15], Min[h, 15]}, PlotStyle -> Opacity[0.5], Filling -> Bottom, 
                              PlotRange -> All, ImageSize -> {300, 300}, SphericalRegion -> True, 
                              ImagePadding -> 15, 
                              PlotLabel -> 
                                  Row[{
                                      Style["Dirichlet Energy: ", Italic,  FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"], 
                                      Style[f,  FontSize -> 12]
                                  }],
                              RegionFunction -> Function[{y, x}, 
                                  (1 <= x <= rX && 1 <= y <= l) (*inner boundary above*)
                               || (rX + rXL <= x <= h && 1 <= y <= l) (*inner boundary below*)
                               || ((rX <= x <= rX + rXL) && 1 <= y <= rY) (*inner boundary left*)
                               || ((rX <= x <= rX + rXL) && rY + rYL + 1 <= y <= l) (*inner boundary right*)], 
                PerformanceGoal -> "Speed"]
            }, Alignment -> Center]
        ]
    ]

(*******************************************************************)
(************************* manipulate ******************************)
(*******************************************************************)
Manipulate[
    SeedRandom[seed];
    matrix  = Table[0, {x, 1 , matrixHeight}, {y, 1, matrixwidth}];

    (*********************** initialize matrix ***********************)
    Switch[boundaryType ,
        "Random[0, d]",
            matrix = randomInitilization[matrix, randomRange, innerBoundaryRandomRange,
                                        innerBoundary, innerBoundaryXIndex, innerBoundaryYIndex,
                                        innerBoundaryXwidth, innerBoundaryYwidth, randomRangeDomain],
        "c + sin(ax + by)",
            matrix = sinInitialization[matrix, innerBoundary, innerBoundaryXIndex, innerBoundaryYIndex, 
                                      innerBoundaryXwidth, innerBoundaryYwidth, aSin, bSin, cSin, 
                                      innerBoundaryA, innerBoundaryB, innerBoundaryC, randomRangeSinDomain],

        TraditionalForm[c + (ax^2 + bx^2)],
            matrix = quadraticInitialization[matrix, innerBoundary, innerBoundaryXIndex, innerBoundaryYIndex, 
                                            innerBoundaryXwidth, innerBoundaryYwidth, quadraA, quadraB, quadraC, 
                                            quaInnerBoundaryA, quaInnerBoundaryB, quaInnerBoundaryC, 
                                            randomRangeQuadraDomain]
    ];

    If[innerBoundary,
    (**************** apply relaxation process 20 times ****************)

        matrixList = 
            NestList[
                relaxationProcessInnerBoundary[#, innerBoundaryXIndex, innerBoundaryYIndex, 
                                                  innerBoundaryXwidth, innerBoundaryYwidth]&, matrix, 20];
        ,
        matrixList = NestList[relaxationProcess[#] &, matrix, 20];
    ];

    matrixList = initializeCorner[matrixList, currentIteration + 1];
    energy = dirichletEnergy[matrixList[[currentIteration + 1]]];

    If[innerBoundary,
        plotGraphInnerBoundary[matrixList[[currentIteration + 1]], energy, 
                              innerBoundaryXIndex, innerBoundaryXwidth, innerBoundaryYIndex, 
                              innerBoundaryYwidth, boundaryType, randomRange, 
                              innerBoundaryRandomRange, aSin, bSin, cSin, innerBoundaryA, 
                              innerBoundaryB, innerBoundaryC, quadraA, quadraB, quadraC, 
                              quaInnerBoundaryA, quaInnerBoundaryB, quaInnerBoundaryC, 
                              currentIteration]
        ,
        plotGraph[matrixList[[currentIteration + 1]], energy, boundaryType, 
                  randomRange, aSin, bSin, cSin, quadraA, quadraB, quadraC, 
                  currentIteration]
    ],
 
(*******************************************************************)
(************************ control area *****************************)
(*******************************************************************)
 
    PaneSelector[{
        "Random[0, d]" -> Column[{
                Control@{{boundaryType, "Random[0, d]",  Style["boundary function", FontSize -> 12, 
                          FontFamily -> "Arial Hebrew Scholar" ]}, 
                          {"Random[0, d]", "c + sin(ax + by)", TraditionalForm[c + (ax^2 + bx^2)]}, 
                          ControlType -> PopupMenu},
                Row[{
                    Spacer[30], 
                    Style["Initialization:   Random[0,  d]", Bold, FontSize -> 12 ]
                }],
                Control@{{randomRangeDomain, 10, "d"}, 1, 50, 1, Appearance -> "Labeled"},
                Control@{{currentIteration, 0, "\!\(\*StyleBox[\"iterations\",\nFontSlant->\"Italic\"]\)"}, 
                        0, 19, 1, AnimationRate -> 2, Appearance -> "Labeled"},
                Button[
                    Style["new randomization", FontSize -> 11, FontFamily -> "Arial Hebrew Scholar" ], 
                                              seed += 1, ImageSize -> {130, 25}, Alignment -> Center],

                (*============ outer boundary ==============*)
                Row[{
                    Spacer[30], Style["outer boundary:   Subscript[c, 0]+(Subscript[a, 0] x^2+Subscript[b, 0] y^2)", FontSize -> 12, Bold]
                }],
                Control@{{randomRange, 25, "Subscript[a, 0]"} , 1, 50, 1, Appearance -> "Labeled", Enabled -> innerBoundary},
                Style["dimensions of outer matrix", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar" ],
                Row[{
                    Spacer[20], 
                    Control@{{matrixHeight, 13, "height"}, Range[8, 50], ControlType -> PopupMenu}, Spacer[20], 
                    Control@ {{matrixwidth, 13, "width"}, Range[6, 50], ControlType -> PopupMenu}
                }],

                (*============ inner boundary ==============*)
                Row[{
                    Control@{{innerBoundary, False, ""}, {False, True}},  
                    Style[" inner boundary:   Random[0,  Subscript[a, 1]]", Bold, FontSize -> 12, FontFamily -> "Arial Hebrew Scholar" ]
                }],
                Control@{{innerBoundaryRandomRange, 50, "Subscript[a, 1]"}, 1, 50, 1, Appearance -> "Labeled", Enabled -> innerBoundary},
                Style["position of inner matrix", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar" ],
                Row[{
                    Spacer[40], Control@ {{innerBoundaryXIndex, 5, "x"}, Range[2, matrixHeight - 1], ControlType -> PopupMenu, 
                    Enabled -> innerBoundary}, Spacer[45], 
                    Control@ {{innerBoundaryYIndex, 5, "y"}, Range[2, matrixwidth - 1], ControlType -> PopupMenu, Enabled -> innerBoundary}
                }],
                Style["dimensions of inner matrix", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar" ],
                Row[{
                    Spacer[20], 
                    Control@ {{innerBoundaryXwidth, 6 , "height"}, Range[0, matrixHeight - 1 - innerBoundaryXIndex], 
                                                                  ControlType -> PopupMenu, Enabled -> innerBoundary}, 
                    Spacer[25],
                    Control@ {{innerBoundaryYwidth, 6 , "width"}, Range[0, matrixwidth - 1 - innerBoundaryYIndex], 
                                                                  ControlType -> PopupMenu, Enabled -> innerBoundary}
                }]
                }, Dividers -> {None, {6 -> RGBColor[0.7764, 0.7725, 0.8], 10 -> RGBColor[0.7764, 0.7725, 0.8]}
            }],

        "c + sin(ax + by)" -> Column[{
                Control@{{boundaryType, "Random[0, d]",  Style["boundary function", FontSize -> 12, 
                          FontFamily -> "Arial Hebrew Scholar" ]}, 
                          {"Random[0, d]", "c + sin(ax + by)", TraditionalForm[c + (ax^2 + bx^2)]}, 
                          ControlType -> PopupMenu},
                Row[{
                    Spacer[30], 
                    Style["Initialization:   Random[0,  d]", Bold, FontSize -> 12 ]
                }],
                Control@{{randomRangeSinDomain, 0, "d"}, 0, 5, 0.1, Appearance -> "Labeled"},
                Control@{{currentIteration, 0, "\!\(\*StyleBox[\"iterations\",\nFontSlant->\"Italic\"]\)"}, 
                          0, 19, 1, AnimationRate -> 2, Appearance -> "Labeled"},
                
                (*============ outer boundary ==============*)
                Row[{
                    Spacer[30], 
                    Style["outer boundary:   Subscript[c, 0]+(Subscript[a, 0] x^2+Subscript[b, 0] y^2)", FontSize -> 12, Bold]
                }],
                Control@{{aSin, .3, "Subscript[a, 0]"}, Range[-4, 4, 0.1], ControlType -> Slider, Appearance -> "Labeled"},
                Control@{{bSin, .4, "Subscript[b, 0]"}, Range[0, 4, 0.1], ControlType -> Slider, Appearance -> "Labeled"},
                Control@{{cSin, 0, "Subscript[c, 0]"}, Range[-4, 4, 0.1], ControlType -> Slider, Appearance -> "Labeled", 
                                                      Enabled -> innerBoundary},
                Style["dimensions of outer matrix", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"],
                Row[{
                    Spacer[20], 
                    Control@{{matrixHeight, 13, "height"}, Range[8, 50], ControlType -> PopupMenu},
                    Spacer[20], 
                    Control@ {{matrixwidth, 13, "width"}, Range[6, 50], ControlType -> PopupMenu}
                }],
                
                (*============ inner boundary ==============*)
                Row[{
                    Control@{{innerBoundary, False, ""}, {False, True}},  
                    Style[" inner boundary:   Subscript[c, 1]+sin(Subscript[a, 1] x + Subscript[b, 1] y)", Bold, FontSize -> 12]
                }],
                Control@{{innerBoundaryA, .2, "Subscript[a, 1]" }, Range[-4, 4, 0.1], ControlType -> Slider, Appearance -> "Labeled", 
                                                                  Enabled -> innerBoundary},
                Control@{{innerBoundaryB, .4, "Subscript[b, 1]"}, Range[0, 4, 0.1], ControlType -> Slider, Appearance -> "Labeled", 
                                                                  Enabled -> innerBoundary},
                Control@{{innerBoundaryC, 0, "Subscript[c, 1]"}, Range[-4, 4, 0.1], ControlType -> Slider, Appearance -> "Labeled", 
                                                                  Enabled -> innerBoundary},
                Style["position of inner matrix", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"],
                Row[{
                    Spacer[40], 
                    Control@ {{innerBoundaryXIndex, 5, "x"}, Range[2, matrixHeight - 1], ControlType -> PopupMenu, Enabled -> innerBoundary}, 
                    Spacer[45], 
                    Control@ {{innerBoundaryYIndex, 5, "y"}, Range[2, matrixwidth - 1], ControlType -> PopupMenu, Enabled -> innerBoundary}
                }],
                Style["dimensions of inner matrix", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"],
                Row[{
                    Spacer[20], 
                    Control@ {{innerBoundaryXwidth, 6 , "height"}, Range[0, matrixHeight - 1 - innerBoundaryXIndex], 
                                                                  ControlType -> PopupMenu, Enabled -> innerBoundary}, 
                    Spacer[25],
                    Control@ {{innerBoundaryYwidth, 6 , "width"}, Range[0, matrixwidth - 1 - innerBoundaryYIndex], 
                                                                  ControlType -> PopupMenu, Enabled -> innerBoundary}
                }]
                }, Dividers -> {None, {5 -> RGBColor[0.7764, 0.7725, 0.8], 11 -> RGBColor[0.7764, 0.7725, 0.8]}
            }],
        TraditionalForm[c + (ax^2 + bx^2)] -> Column[{
            Control@{{boundaryType, "Random[0, d]",  
            Style["boundary function", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"]}, 
            {"Random[0, d]", "c + sin(ax + by)", TraditionalForm[c + (ax^2 + bx^2)]}, 
            ControlType -> PopupMenu},
            Row[{
                Spacer[30], 
                Style["Initialization:   Random[0,  d]", Bold, FontSize -> 12 ]}],
                Control@{{randomRangeQuadraDomain, 0, "d"}, 0, 500, 10, Appearance -> "Labeled"},
                Control@{{currentIteration, 0, "\!\(\*StyleBox[\"iterations\",\nFontSlant->\"Italic\"]\)"}, 
                        0, 19, 1, AnimationRate -> 2, Appearance -> "Labeled"},
            
            (*============ outer boundary ==============*)
            Row[{
                Spacer[30], 
                Style["outer boundary:   Subscript[c, 0]+(Subscript[a, 0] x^2+Subscript[b, 0] y^2)", FontSize -> 12, Bold]
            }],
            Control@{{quadraA, 4, "Subscript[a, 0]"}, Range[-4, 4, 0.1], ControlType -> Slider, Appearance -> "Labeled"},
            Control@{{quadraB, 0, "Subscript[b, 0]"}, Range[0, 4, 0.1], ControlType -> Slider, Appearance -> "Labeled"},
            Control@{{quadraC, .4, "Subscript[c, 0]"}, Range[-4, 4, 0.1], ControlType -> Slider, Appearance -> "Labeled", 
                                                      Enabled -> innerBoundary},
            Style["dimensions of outer matrix", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"],
            Row[{
                Spacer[20], 
                Control@{{matrixHeight, 13, "height"}, Range[8, 50], ControlType -> PopupMenu}, 
                Spacer[20], 
                Control@ {{matrixwidth, 13, "width"}, Range[6, 50], ControlType -> PopupMenu}
            }],
            
            (*============ inner boundary ==============*)
            Row[{
                Control@{{innerBoundary, False, ""}, {False, True} },  
                Style[" inner boundary:   Subscript[c, 1]+(Subscript[a, 1] x^2+Subscript[b, 1] y^2)", Bold, FontSize -> 12]
            }],
            Control@{{quaInnerBoundaryA, .2, "Subscript[a, 1]"}, Range[-4, 4, 0.1], ControlType -> Slider, Appearance -> "Labeled", 
                                                                  Enabled -> innerBoundary},
            Control@{{quaInnerBoundaryB, .4, "Subscript[b, 1]"}, Range[0, 4, 0.1], ControlType -> Slider, Appearance -> "Labeled", 
                                                                  Enabled -> innerBoundary},
            Control@{{quaInnerBoundaryC, .4, "Subscript[c, 1]"}, Range[-4, 4, 0.1], ControlType -> Slider, Appearance -> "Labeled",
                                                                  Enabled -> innerBoundary},
            Style["position of inner matrix", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar"],
            Row[{
                Spacer[40], 
                Control@ {{innerBoundaryXIndex, 5, "x"}, Range[2, matrixHeight - 1], ControlType -> PopupMenu, Enabled -> innerBoundary}, 
                Spacer[45], 
                Control@ {{innerBoundaryYIndex, 5, "y"}, Range[2, matrixwidth - 1], ControlType -> PopupMenu, Enabled -> innerBoundary}
            }],
            Style["dimensions of inner matrix", FontSize -> 12, FontFamily -> "Arial Hebrew Scholar" ],
            Row[{
                Spacer[20], 
                Control@ {{innerBoundaryXwidth, 6 , "height"}, Range[0, matrixHeight - 1 - innerBoundaryXIndex], 
                                                              ControlType -> PopupMenu, Enabled -> innerBoundary}, 
                Spacer[25],
                Control@ {{innerBoundaryYwidth, 6 , "width"}, Range[0, matrixwidth - 1 - innerBoundaryYIndex], 
                                                              ControlType -> PopupMenu, Enabled -> innerBoundary}}]
            }, Dividers -> {None, {5 -> RGBColor[0.7764, 0.7725, 0.8], 11 -> RGBColor[0.7764, 0.7725, 0.8]}
        }]
    }, Dynamic@boundaryType],
    {seed, 1, Range[0, 100], ControlType -> None},
    ControlPlacement -> Left,
    TrackedSymbols :> True,
    SaveDefinitions -> True,
    SynchronousInitialization -> False
]
