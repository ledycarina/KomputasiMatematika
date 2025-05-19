function genes = create_genes (led)
random_number = randi ([32,126],1,led);
genes = char(random_number);
end