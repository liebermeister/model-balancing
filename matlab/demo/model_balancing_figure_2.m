% Generate Figure 2 from "Model balancing: a search for in-vivo kinetic constants and consistent metabolic states", Liebermeister and Noor

v = 1;
s = 1;
p = 0.1;
keq1 = 1;
keq2 = 1;
kcat1 = 5;
kcat2 = 1;

a = 10.^[-0.99:0.01:-0.01];
theta1 = log(keq1) - log(a/s);
theta2 = log(keq2) - log(p./a);

x  = log(a);
y1 = log( v./(kcat1 * (1-exp(-theta1))));
y2 = log( v./(kcat1 * (1-exp(-theta2))));

y1_data  = 0.2;
y1_sigma = 0.2;

y2_data  = -0.3;
y2_sigma = 0.1;

y1_score = 0.5 * (y1-y1_data).^2/y1_sigma^2;
y2_score = 0.5 * (y2-y2_data).^2/y2_sigma^2;

y1_score_convex = y1_score .* double(y1>y1_data);
y2_score_convex = y2_score .* double(y2>y2_data);
alpha = 0.2;
y1_score_alpha  = 0.2 * y1_score + 0.8 * y1_score_convex;
y2_score_alpha  = 0.2 * y2_score + 0.8 * y2_score_convex;

rest = 20*(x-(-1.2)).^2;

figure(1); clf;
subplot(4,1,1); plot(x,y1,'r'); axis([min(x), max(x), -2,1]); hold on;
plot([min(x),max(x)],y1_data*[1,1],'r-');
plot([min(x),max(x)],[y1_data+y1_sigma]*[1,1],'r--');
plot([min(x),max(x)],[y1_data-y1_sigma]*[1,1],'r--');
ylabel('ln e = y(x)')
subplot(4,1,2); plot(x,y1_score,'k');axis([min(x), max(x), 0,100]); legend('alpha=1');
ylabel('Py(x)')
subplot(4,1,3); plot(x,y1_score_alpha,'b');axis([min(x), max(x), 0,100]); legend('alpha=0.2');
ylabel('Py(x)')
subplot(4,1,4); plot(x,y1_score_convex,'c'); axis([min(x), max(x), 0,100]); legend('alpha=0');
ylabel('Py(x)')
xlabel('x = ln c');

figure(2); clf;
subplot(4,1,1); plot(x,y2,'r.'); axis([min(x), max(x), -2,1]);  hold on;
plot([min(x),max(x)],y2_data*[1,1],'r-');
plot([min(x),max(x)],[y2_data+y2_sigma]*[1,1],'r--');
plot([min(x),max(x)],[y2_data-y2_sigma]*[1,1],'r--');
ylabel('ln e = y(x)')
subplot(4,1,2); plot(x,y2_score,'k.');axis([min(x), max(x), 0,100]); legend('alpha=1');
ylabel('Py(x)')
subplot(4,1,3); plot(x,y2_score_alpha,'b.');axis([min(x), max(x), 0,100]); legend('alpha=0.2');
ylabel('Py(x)')
subplot(4,1,4); plot(x,y2_score_convex,'c.'); axis([min(x), max(x), 0,100]); legend('alpha=0');
ylabel('Py(x)')
xlabel('x = ln c');

figure(3); clf;

subplot(3,1,1); 
plot(x,y1_score,'k-');hold on; 
plot(x,y2_score,'k.');        
plot(x,y1_score+y2_score,'k','Linewidth',2);        
plot(x,y1_score+y2_score+rest,'k--','Linewidth',2);        
axis([min(x), max(x), 0,100]); legend('Py,alpha=1','Total score');
ylabel('Py(x)')

subplot(3,1,2); 
plot(x,y1_score_alpha,'b-');hold on; 
plot(x,y2_score_alpha,'b.');  
plot(x,y1_score_alpha+y2_score_alpha,'b','Linewidth',2);  
plot(x,y1_score_alpha+y2_score_alpha+rest,'b--','Linewidth',2);  
axis([min(x), max(x), 0,100]); legend('Py,alpha=0.2','Total score');
ylabel('Py(x)')

subplot(3,1,3); 
plot(x,y1_score_convex,'c-');hold on; 
plot(x,y2_score_convex,'c.'); 
plot(x,y1_score_convex+y2_score_convex,'c','Linewidth',2); 
plot(x,y1_score_convex+y2_score_convex+rest,'c--','Linewidth',2); 
axis([min(x), max(x), 0,100]); legend('Py,alpha=1','Total score');
ylabel('Py(x)')
xlabel('x = ln c');
