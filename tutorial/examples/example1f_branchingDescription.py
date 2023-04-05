

'''
In this new branch, the definition of organ lateral is more flexible. 
It is possible to define for each node the number, organ type, and subtype of the laterals
TODO add tutorial

for backward compatibility:
	- successors of stems, with subtype 2 and no specified organ type will be set as leaves
	- default value of "organType" of the successor will be identical to the organType of the parent organ
	- default value for probability is 1
	- default value of number of lateral per successor rule is 1

TODO: add example
'''