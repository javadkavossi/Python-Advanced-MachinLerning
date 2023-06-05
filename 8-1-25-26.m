% تمرین ۸.۲۵ - الگوریتم شیب‌دارترین با روش secant

% تابع هدف (مثال ۸.۱)
function fval = f(x)
    % تابع هدف
    % ورودی:
    % x: بردار ورودی
    % خروجی:
    % fval: مقدار تابع هدف برای بردار ورودی
    
    % محاسبه مقدار تابع هدف
    fval = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
end

% الگوریتم شیب‌دارترین با روش secant
function x_star = gradientDescentSecant(x0, epsilon)
    % الگوریتم شیب‌دارترین با روش secant
    % ورودی:
    % x0: بردار شرط اولیه
    % epsilon: مقدار کوچک جهت اتمام الگوریتم
    % خروجی:
    % x_star: بردار نقطه پایانی
    
    % تعریف متغیرهای اولیه
    x_prev = x0;
    f_prev = f(x_prev);
    
    % محاسبه گرادیان
    grad_f_prev = [400*(x_prev(1)^3 - x_prev(1)*x_prev(2)) + 2*x_prev(1) - 2; 200*(x_prev(2) - x_prev(1)^2)];
    
    % محاسبه متغیرهای بعدی
    x = x_prev - grad_f_prev * epsilon;
    f_curr = f(x);
    
    % حلقه اصلی الگوریتم
    while norm(grad_f_prev) > 1e-6
        x_temp = x;
        f_temp = f_curr;
        
        % محاسبه گرادیان
        grad_f_temp = [400*(x_temp(1)^3 - x_temp(1)*x_temp(2)) + 2*x_temp(1) - 2; 200*(x_temp(2) - x_temp(1)^2)];
        
        % محاسبه متغیر بعدی با روش secant
        x = x_temp - (grad_f_temp * (x_temp - x_prev)) / (grad_f_temp - grad_f_prev) * (x_temp - x_prev);
        
        % به‌روزرسانی متغیرها
        x_prev = x_temp;
        f_prev = f_temp;
        grad_f_prev = grad_f_temp;
        f_curr = f(x);
    end
    
    % تخمین نهایی
    x_star = x;
end

% تمرین ۸.۲۵ - آزمایش برنامه با شرط اولیه [-4, 5, 1]
x0 = [-4; 5; 1]; % شرط اولیه
epsilon = 1e-6; % مقدار کوچک جهت اتمام الگوریتم
x_star = gradientDescentSecant(x0, epsilon);
disp('نقطه پایانی:');
disp(x_star);
disp('مقدار تابع هدف در نقطه پایانی:');
disp(f(x_star));

% تمرین ۸.۲۶ - تابع هدف جدید و آزمایش برنامه با شرط اولیه [0, 0]
function fval = f_new(x)
    % تابع هدف جدید
    % ورودی:
    % x: بردار ورودی
    % خروجی:
    % fval: مقدار تابع هدف جدید برای بردار ورودی
    
    % محاسبه مقدار تابع هدف جدید
    fval = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
end

x0_new = [0; 0]; % شرط اولیه
x_star_new = gradientDescentSecant(x0_new, epsilon);
disp('نقطه پایانی (تابع هدف جدید):');
disp(x_star_new);
disp('مقدار تابع هدف در نقطه پایانی (تابع هدف جدید):');
disp(f_new(x_star_new));
