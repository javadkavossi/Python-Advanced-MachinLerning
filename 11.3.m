% تمرین 11.3 - الگوریتم DFP برای مینیمایز کردن تابع f(a)

% تابع هدف
function fval = f(a)
    % تابع هدف
    % ورودی:
    % a: بردار ورودی
    % خروجی:
    % fval: مقدار تابع هدف برای بردار ورودی
    
    % محاسبه مقدار تابع هدف
    fval = 20 + (a(1) + 1)^2 + (a(2) - 1)^2 - 10*(cos(2*pi*(a(1) + 1)) + cos(2*pi*(a(2) - 1)));
end

% الگوریتم DFP
function a_star = minimizeDFP(a0, epsilon)
    % الگوریتم DFP برای مینیمایز کردن تابع f(a)
    % ورودی:
    % a0: بردار شرط اولیه
    % epsilon: مقدار کوچک جهت اتمام الگوریتم
    % خروجی:
    % a_star: بردار نقطه پایانی
    
    % تعریف متغیرهای اولیه
    a_prev = a0;
    f_prev = f(a_prev);
    
    % محاسبه گرادیان
    grad_f_prev = [2*(a_prev(1) + 1) + 20*pi*sin(2*pi*(a_prev(1) + 1)); 2*(a_prev(2) - 1) - 20*pi*sin(2*pi*(a_prev(2) - 1))];
    
    % ماتریس H0
    H0 = eye(2);
    
    % حلقه اصلی الگوریتم
    while norm(grad_f_prev) > epsilon
        % محاسبه جهت جستجو
        p = -H0 * grad_f_prev;
        
        % تابع خط جستجو
        alpha = lineSearch(a_prev, p);
        
        % به‌روزرسانی متغیرها
        a = a_prev + alpha * p;
        f_curr = f(a);
        
        % محاسبه گرادیان
        grad_f_curr = [2*(a(1) + 1) + 20*pi*sin(2*pi*(a(1) + 1)); 2*(a(2) - 1) - 20*pi*sin(2*pi*(a(2) - 1))];
        
        % محاسبه ماتریس H
        delta_a = a - a_prev;
        delta_grad_f = grad_f_curr - grad_f_prev;
        H = H0 + (delta_a * delta_a') / (delta_a' * delta_grad_f) - (H0 * delta_grad_f * delta_grad_f' * H0) / (delta_grad_f' * H0 * delta_grad_f);
        
        % به‌روزرسانی متغیرهای قبلی
        a_prev = a;
        f_prev = f_curr;
        grad_f_prev = grad_f_curr;
        H0 = H;
    end
    
    % تخمین نهایی
    a_star = a_prev;
end

% تابع جستجوی خط
function alpha = lineSearch(a, p)
    % تابع جستجوی خط با استفاده از روش secant
    % ورودی:
    % a: بردار شرط اولیه
    % p: جهت جستجو
    % خروجی:
    % alpha: مقدار پارامتر جستجو
    
    alpha0 = 0;
    alpha1 = 0.1;
    epsilon = 1e-6;
    
    while abs(f(a + alpha1*p) - f(a + alpha0*p)) > epsilon
        alpha2 = alpha1 - (alpha1 - alpha0) * f(a + alpha1*p) / (f(a + alpha1*p) - f(a + alpha0*p));
        alpha0 = alpha1;
        alpha1 = alpha2;
    end
    
    alpha = alpha1;
end

% اجرای الگوریتم با شرط اولیه (0) = [0, 0]
a0 = [0; 0];
epsilon = 1e-6;
a_star = minimizeDFP(a0, epsilon);
disp('نقطه پایانی:');
disp(a_star);
disp('مقدار تابع هدف در نقطه پایانی:');
disp(f(a_star));
