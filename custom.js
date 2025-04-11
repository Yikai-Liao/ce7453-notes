/**
 * 在 mdBook 页面导航栏插入“下载 PDF”按钮
 */
(function () {
    function insertPdfButton() {
        // 检查是否已插入，避免重复
        if (document.getElementById('download-pdf-btn')) return;

        // 创建按钮
        var btn = document.createElement('a');
        btn.id = 'download-pdf-btn';
        btn.href = '/CE7453.pdf';
        btn.target = '_blank';
        btn.rel = 'noopener';
        btn.textContent = '下载 PDF';
        btn.style.marginLeft = '16px';
        btn.style.padding = '6px 12px';
        btn.style.background = '#007bff';
        btn.style.color = '#fff';
        btn.style.borderRadius = '4px';
        btn.style.textDecoration = 'none';
        btn.style.fontSize = '14px';

        // 定位到导航栏右侧
        var nav = document.querySelector('.nav-links');
        if (nav) {
            nav.appendChild(btn);
        } else {
            // fallback: 插入到 body
            document.body.appendChild(btn);
        }
    }

    // mdBook 页面加载后执行
    if (window.addEventListener) {
        window.addEventListener('DOMContentLoaded', insertPdfButton, false);
    } else if (window.attachEvent) {
        window.attachEvent('onload', insertPdfButton);
    }
})();