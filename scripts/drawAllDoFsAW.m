format = '.jpg';
res = 10;
h = 1e-10;
% h = 1;

for j = 1:24	
% 	for k = 1:2
k = 1;
		clf
		drawDoF(j,k,res,h);
		axis off;
		grid off;
		k = 4;
		filename = ['DoF_' num2str(j) '_Comp_' num2str(k)];
		saveas(gcf,[filename format]);
% 	end
end

% for j = 1:24	
% % 	j = 9;
% 	k = 1;
% % 	for k = 1:3
% 		clf
% 		drawDoFnormal(j,k,res);
% 		axis off;
% 		grid off;
% 		hold off;
% 		colorbar
% 		
% 		k = 100;
% 		filename = ['DoF_' num2str(j) '_Normal_' num2str(k)];
% 		title(['sigma-' num2str(j)  '-Normal-' num2str(k)] );
% 		saveas(gcf,[filename format]);
% % 	end
% end
