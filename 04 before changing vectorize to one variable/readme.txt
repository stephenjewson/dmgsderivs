1
Everything up to 52 works, and gives correct ppm using test8.
So that's a good start.
Note that frechet_k1 was a pain, because of the parameter ordering switching round, which is my stupid fault.

2
I've now added f1f and f2f to non-predictor models, but haven't tested it.
And I've half added p1f and p2f (just the deriv, not the wrapper), but haven't tested it yet.
Next steps would be:

-test f1f and f2f in some real ppm cases

-if they work, extend them to predictor models


-finish p1f and p2f for non-predictors

-test them

-if they work, extend them to predictor models

At that point, I'm done with aderivs, I think.

s 