function result = my_exp(base,exponent)
	if (exponent == 2) 
		result = base*base;
        elseif (exponent == 0.5) 
		result = sqrt(base);
	else 
		error('This is weird');
	end
end
