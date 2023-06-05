% تمرین 11.17 - الگوریتم شبه نیوتن برای توابع عمومی با روش secant

% تابع هدف (مثال 9.4)
function fval = f(x)
    % تابع هدف
    % ورودی:
    % x: بردار ورودی
    % خروجی:
    % fval: مقدار تابع هدف برای بردار ورودی
    
    % محاسبه مقدار تابع هدف
    fval = (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;
end

% الگوریتم شبه نیوتن
function x_star = quasiNewton(f, x0)
    % الگوریتم شبه نیوتن برای توابع عمومی
    % ورودی:
    % f: تابع هدف
    % x0: بردار شرط اولیه
    % خروجی:
    % x_star: بردار نقطه پایانی
    
    % تعریف متغیرهای اولیه
    x_prev = x0;
    f_prev = f(x_prev);
    grad_f_prev = gradient(f, x_prev);
    H_prev = eye(length(x0));
    
    % حلقه اصلی الگوریتم
    while norm(grad_f_prev) > 1e-6
        % جهت به‌روزرسانی
        if mod(iter, 6) == 0
            p = -grad_f_prev;
        else
            p = -H_prev * grad_f_prev;
        end
        
        % تابع خط جستجو
        alpha = lineSearch(f, x_prev, p);
        
        % به‌روزرسانی متغیرها
        x = x_prev + alpha * p;
        f_curr = f(x);
        grad_f_curr = gradient(f, x);
        
        % به‌روزرسانی ماتریس H
        s = x - x_prev;
        y = grad_f_curr - grad_f_prev;
        rho = 1 / (y * s);
        H = (eye(length(x0)) - rho * s * y) * H_prev * (eye(length(x0)) - rho * y * s') + rho * s * s';
        
        % به‌روزرسانی متغیرهای قبلی
        x_prev = x;
        f_prev = f_curr;
        grad_f_prev = grad_f_curr;
        H_prev = H;
    end
    
    % تخمین نهایی
    x_star = x_prev;
end

% تابع جستجوی خط
function alpha = lineSearch(f, x, p)
    % تابع جستجوی خط با استفاده از روش secant
    % ورودی:
    % f: تابع هدف
    % x: بردار شرط اولیه
    % p: جهت جستجو
    % خروجی:
    % alpha: مقدار پارامتر جستجو
    
    alpha0 = 0;
    alpha1 = 0.1;
    epsilon = 1e-6;
    
    while abs(f(x + alpha1 * p) - f(x + alpha0 * p)) > epsilon
        alpha2 = alpha1 - (alpha1 - alpha0) * f(x + alpha1 * p) / (f(x + alpha1 * p) - f(x + alpha0 * p));
        alpha0 = alpha1;
        alpha1 = alpha2;
    end
    
    alpha = alpha1;
end

% اجرای الگوریتم با شرط اولیه (0) = [-2،2]
x0 = [-2; 2];
x_star = quasiNewton(@f, x0);
disp('نقطه پایانی:');
disp(x_star);
disp('مقدار تابع هدف در نقطه پایانی:');
disp(f(x_star));
