# FurnaceTest

White furnace test and weak white furnace test implementation for Cook-Torrance with GGX distribution and Smith (height-correlated) geometry term and square roughness remap. 

White furnace test is basically a lighting integral against a constant white environment. It should always evaluate to 1 if there is no energy loss or gain. In this case we can see substantial amount of energy loss, especially for high roughness values.

More info:
Eric Heitz - "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs", JCGT 2014
http://jcgt.org/published/0003/02/03/